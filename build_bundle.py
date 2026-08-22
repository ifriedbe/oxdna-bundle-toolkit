"""Build an N-helix DNA duplex bundle for oxDNA, as an idealized proxy
test of whether bundling (the standard DNA-nanotech fix for a floppy
single duplex) closes the stiffness gap between a single duplex and a
usable rigid nanomechanical element. In this project's own testing, a
single 20bp duplex measured ~3,100x too soft for a target stiffness of
0.41 N/m; see the project README for the full story.

IMPORTANT HONESTY NOTE: this does NOT build real covalent crossovers
(the actual mechanism DNA origami uses to rigidly lock adjacent helices
together -- a strand physically crosses from one helix to the next at
specific base-paired, torsionally-registered points). Building real
crossovers correctly requires helical-register-aware placement that a
dedicated origami CAD tool (cadnano, etc.) handles; hand-building it
here risked an unphysical, silently-broken structure. Instead, this
script uses oxDNA's `mutual_trap` external force to hold the six
helices' tip ends in a rigid hexagonal arrangement relative to each
other -- an IDEALIZED UPPER BOUND on what real crosslinking could
achieve (traps enforce perfect rigidity at the linked points with zero
extra flexibility, which real crossovers do not). If even this
best-case proxy fails to reach the needed stiffness, that's strong
evidence bundling alone can't close the gap. If it succeeds, it's a
necessary-but-not-sufficient positive signal -- a real crossover-based
structure would likely be less rigid than this idealized case, not
more.

Six 20bp duplexes (same sequence as the single-duplex test, non-
palindromic so copies don't cross-hybridize), arranged in a regular
hexagon, circumradius 2.5nm (~2.935 oxDNA length units), all parallel
along +z -- a standard DNA-nanotech helix-helix spacing.
"""

import os
import numpy as np

OXDNA_LENGTH_NM = 0.8518

# ---- read physical constants from oxDNA's model.h, same as generate-sa.py ----
def import_model_constants():
    consts = {}
    model_path = os.path.join(
        os.path.dirname(__file__), "oxdna_src", "src", "model.h"
    )
    with open(model_path) as f:
        for line in f:
            line = line.split("//")[0].strip()
            if not line.startswith("#define "):
                continue
            rest = line[len("#define "):].strip().split(" ", 1)
            if len(rest) != 2:
                continue
            key, val = [x.strip() for x in rest]
            val = val.replace("f", "").replace("PI", str(np.pi))
            try:
                consts[key] = eval(val, {}, consts)
            except Exception:
                pass
    return consts

_C = import_model_constants()
POS_BASE = _C["POS_BASE"]
CM_CENTER_DS = POS_BASE + 0.2
BASE_BASE = 0.3897628551303122

number_to_base = {0: "A", 1: "G", 2: "C", 3: "T"}
base_to_number = {"A": 0, "G": 1, "C": 2, "T": 3}


def get_rotation_matrix(axis, angle):
    axis = np.array(axis, dtype=float)
    axis /= np.sqrt(np.dot(axis, axis))
    ct, st, olc = np.cos(angle), np.sin(angle), 1.0 - np.cos(angle)
    x, y, z = axis
    return np.array([
        [olc*x*x+ct, olc*x*y-st*z, olc*x*z+st*y],
        [olc*x*y+st*z, olc*y*y+ct, olc*y*z-st*x],
        [olc*x*z-st*y, olc*y*z+st*x, olc*z*z+ct],
    ])


def generate_duplex(bp, sequence, start_pos, dir_vec, perp):
    positions, a1s, a3s = [], [], []
    dir_vec = np.array(dir_vec, dtype=float)
    dir_vec /= np.linalg.norm(dir_vec)
    v1 = np.array(perp, dtype=float)
    v1 -= dir_vec * np.dot(dir_vec, v1)
    v1 /= np.linalg.norm(v1)

    R = get_rotation_matrix(dir_vec, np.pi / 180.0 * 35.9)  # ~1bp twist

    a1, a3, rb = v1, dir_vec, np.array(start_pos, dtype=float)
    for i in range(bp):
        rcdm = rb - CM_CENTER_DS * a1
        positions.append(rcdm); a1s.append(a1); a3s.append(a3)
        if i != bp - 1:
            a1 = R @ a1
            rb = rb + a3 * BASE_BASE

    # complementary strand, antiparallel
    a1, a3 = -a1, -dir_vec
    Rt = R.T
    for _ in range(bp):
        rcdm = rb - CM_CENTER_DS * a1
        positions.append(rcdm); a1s.append(a1); a3s.append(a3)
        a1 = Rt @ a1
        rb = rb + a3 * BASE_BASE

    return positions, a1s, a3s


def main():
    import sys
    bp = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    n_helix = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    # number of EXTRA crosslink rings strictly between root and tip,
    # evenly spaced -- approximates real DNA-origami crossovers repeated
    # periodically along the length, instead of only at the two ends.
    # 0 (default) reproduces the earlier tip-only-crosslink behavior.
    n_extra_rings = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    # trap/crosslink stiffness in oxDNA reduced units -- diagnostic knob
    # to test whether finite trap compliance itself (not bundle geometry)
    # is capping the measured tip rigidity.
    trap_stiff = float(sys.argv[4]) if len(sys.argv) > 4 else 5.0
    seq_letters = "GACTGATCGATCGGATCTGA"[:bp]
    seq = [base_to_number[c] for c in seq_letters]

    # keep nearest-neighbor helix spacing at a realistic 2.5nm regardless
    # of ring size -- more helices form a WIDER ring, not a more cramped
    # one, which is the geometrically-motivated way to gain bending
    # stiffness (~r^2), similar to why a hollow tube resists bending
    # better than a solid rod of the same material amount packed tight.
    neighbor_spacing_nm = 2.5
    circumradius_nm = neighbor_spacing_nm / (2 * np.sin(np.pi / n_helix))
    circumradius = circumradius_nm / OXDNA_LENGTH_NM

    box_side = max(30.0, 2 * circumradius + 15.0)
    center = np.array([box_side / 2, box_side / 2, 5.0])
    dir_vec = np.array([0.0, 0.0, 1.0])

    all_positions, all_a1s, all_a3s = [], [], []
    helix_root_idx = []   # (top_start, bottom_end) global indices, per helix
    helix_tip_idx = []    # (top_end, bottom_start) global indices, per helix
    # ring_fractions: interior crosslink rings (excludes root=0.0, tip=1.0,
    # which are handled separately as before), evenly spaced in (0, 1).
    ring_fractions = [(i + 1) / (n_extra_rings + 1) for i in range(n_extra_rings)]
    helix_ring_idx = [[] for _ in ring_fractions]  # per ring: list of (top_idx, bottom_idx) per helix

    for h in range(n_helix):
        angle = 2 * np.pi * h / n_helix
        offset = circumradius * np.array([np.cos(angle), np.sin(angle), 0.0])
        start_pos = center + offset
        perp = offset / np.linalg.norm(offset)  # radially outward, arbitrary valid perp

        pos, a1s, a3s = generate_duplex(bp, seq, start_pos, dir_vec, perp)

        base = len(all_positions)
        all_positions.extend(pos); all_a1s.extend(a1s); all_a3s.extend(a3s)

        top_start = base + 0
        top_end = base + (bp - 1)
        bottom_start = base + bp
        bottom_end = base + (2 * bp - 1)

        helix_root_idx.append((top_start, bottom_end))
        helix_tip_idx.append((top_end, bottom_start))

        for ring_i, frac in enumerate(ring_fractions):
            top_local = min(bp - 1, max(0, round(frac * (bp - 1))))
            bottom_local = bp + (bp - 1 - top_local)
            helix_ring_idx[ring_i].append((base + top_local, base + bottom_local))

    nnucl = len(all_positions)
    nstrands = 2 * n_helix

    # ---- topology file ----
    with open("generated.top", "w") as out:
        print(nnucl, nstrands, file=out)
        nt_count, strand_count = 0, 1
        for h in range(n_helix):
            top_seq = [number_to_base[s] for s in seq]
            bot_seq_nums = [3 - s for s in seq]
            bot_seq_nums.reverse()
            bot_seq = [number_to_base[s] for s in bot_seq_nums]

            for strand in (top_seq, bot_seq):
                print(strand_count, strand[0], -1, nt_count + 1, file=out)
                nt_count += 1
                for i in range(1, len(strand) - 1):
                    print(strand_count, strand[i], nt_count - 1, nt_count + 1, file=out)
                    nt_count += 1
                print(strand_count, strand[-1], nt_count - 1, -1, file=out)
                nt_count += 1
                strand_count += 1

    # ---- configuration file ----
    with open("generated.dat", "w") as out:
        print("t = 0", file=out)
        print("b = ", box_side, box_side, box_side, file=out)
        print("E = 0. 0. 0.", file=out)
        for p, a1, a3 in zip(all_positions, all_a1s, all_a3s):
            print(p[0], p[1], p[2], a1[0], a1[1], a1[2], a3[0], a3[1], a3[2],
                  0., 0., 0., 0., 0., 0., file=out)

    # ---- external forces: anchor all root ends, cross-link all tip ends ----
    with open("ext.dat", "w") as out:
        # anchor every root-end particle (both strands) at its initial position
        for (top_start, bottom_end) in helix_root_idx:
            for idx in (top_start, bottom_end):
                x, y, z = all_positions[idx]
                out.write(f"{{\ntype = trap\nparticle = {idx}\nstiff = {trap_stiff}\n"
                          f"pos0 = {x}, {y}, {z}\nrate = 0.\ndir = 0.,0.,0.\n}}\n\n")

        def crosslink_ring(ring_idx_list, label):
            for h in range(n_helix):
                h_next = (h + 1) % n_helix
                pa1, pa2 = ring_idx_list[h]
                pb1, pb2 = ring_idx_list[h_next]
                for pa, pb in ((pa1, pb1), (pa2, pb2)):
                    ra = np.array(all_positions[pa])
                    rb = np.array(all_positions[pb])
                    r0 = np.linalg.norm(ra - rb)
                    out.write(f"{{\ntype = mutual_trap\nparticle = {pa}\nref_particle = {pb}\n"
                              f"stiff = {trap_stiff}\nr0 = {r0:.6f}\nPBC=0\n}}\n\n")
                    out.write(f"{{\ntype = mutual_trap\nparticle = {pb}\nref_particle = {pa}\n"
                              f"stiff = {trap_stiff}\nr0 = {r0:.6f}\nPBC=0\n}}\n\n")

        # cross-link adjacent helices' tip ends (idealized rigid crosslink proxy)
        crosslink_ring(helix_tip_idx, "tip")

        # cross-link adjacent helices at interior points along the length,
        # approximating periodic real-world crossover spacing
        for ring_i, ring_positions in enumerate(helix_ring_idx):
            crosslink_ring(ring_positions, f"ring{ring_i}")

    # record tip particle indices for the analysis script
    with open("tip_particles.txt", "w") as out:
        for h in range(n_helix):
            out.write(f"{helix_tip_idx[h][0]} {helix_tip_idx[h][1]}\n")

    print(f"Built {n_helix}-helix bundle: {nnucl} nucleotides, {nstrands} strands.")
    print(f"Circumradius: {circumradius_nm} nm ({circumradius:.3f} oxDNA units)")
    print("Wrote generated.top, generated.dat, ext.dat, tip_particles.txt")


if __name__ == "__main__":
    main()
