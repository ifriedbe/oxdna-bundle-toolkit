"""Build a hinge-lever bistable candidate: a mechanical element with
two well-separated states (e.g. "bump" vs "flat"), rather than a stiff
continuously-positioned cantilever. Useful when the actual requirement
is reliable two-state readability (a signal-to-noise separation
criterion) rather than absolute sub-angstrom positional stiffness --
a materially easier target that's often conflated with the stiffer
one. This script reuses the
validated multi-helix bundle geometry from build_bundle.py (real,
measured noise sigma from that work) and adds a biasing external force
that pulls the tip toward an OFFSET target position, representing the
"activated" (tilted/bump) state after a strand-displacement signal
locks the hinge there. Running with and without the bias measures BOTH
states directly: their mean-position separation (the real throw) and
the noise (sigma) within each state -- a direct simulated test of the
bistable-readability design, not just algebra.

Reuses generate_duplex() and the hexagonal-ring bundle construction
from build_bundle.py; only the external-forces section differs (adds
an optional single-particle `trap` pulling the tip toward an offset
target, in addition to the root anchors and tip crosslink ring).
"""

import os
import sys
import numpy as np

OXDNA_LENGTH_NM = 0.8518

sys.path.insert(0, os.path.dirname(__file__))
from build_bundle import (
    import_model_constants, generate_duplex, base_to_number, number_to_base,
)

_C = import_model_constants()


def main():
    bp = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    n_helix = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    bias_nm = float(sys.argv[3]) if len(sys.argv) > 3 else 3.0  # target lateral offset for the "activated" state
    trap_stiff = float(sys.argv[4]) if len(sys.argv) > 4 else 5.0
    # number of tip particles that get the bias lock (realistic strand-
    # displacement lock = 1-2 hybridization points, i.e. one helix's tip
    # pair = 2 particles; default 0 means "all tip particles", matching
    # the earlier, more idealized/optimistic test.
    n_lock_particles = int(sys.argv[5]) if len(sys.argv) > 5 else 0

    seq_letters = "GACTGATCGATCGGATCTGA"[:bp]
    seq = [base_to_number[c] for c in seq_letters]

    neighbor_spacing_nm = 2.5
    circumradius_nm = neighbor_spacing_nm / (2 * np.sin(np.pi / n_helix))
    circumradius = circumradius_nm / OXDNA_LENGTH_NM

    box_side = max(30.0, 2 * circumradius + 15.0)
    center = np.array([box_side / 2, box_side / 2, 5.0])
    dir_vec = np.array([0.0, 0.0, 1.0])

    all_positions, all_a1s, all_a3s = [], [], []
    helix_root_idx = []
    helix_tip_idx = []

    for h in range(n_helix):
        angle = 2 * np.pi * h / n_helix
        offset = circumradius * np.array([np.cos(angle), np.sin(angle), 0.0])
        start_pos = center + offset
        perp = offset / np.linalg.norm(offset)

        pos, a1s, a3s = generate_duplex(bp, seq, start_pos, dir_vec, perp)

        base = len(all_positions)
        all_positions.extend(pos); all_a1s.extend(a1s); all_a3s.extend(a3s)

        top_start = base + 0
        top_end = base + (bp - 1)
        bottom_start = base + bp
        bottom_end = base + (2 * bp - 1)

        helix_root_idx.append((top_start, bottom_end))
        helix_tip_idx.append((top_end, bottom_start))

    nnucl = len(all_positions)
    nstrands = 2 * n_helix

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

    with open("generated.dat", "w") as out:
        print("t = 0", file=out)
        print("b = ", box_side, box_side, box_side, file=out)
        print("E = 0. 0. 0.", file=out)
        for p, a1, a3 in zip(all_positions, all_a1s, all_a3s):
            print(p[0], p[1], p[2], a1[0], a1[1], a1[2], a3[0], a3[1], a3[2],
                  0., 0., 0., 0., 0., 0., file=out)

    # compute tip centroid (undeflected/"flat" state reference) to derive
    # a lateral bias direction perpendicular to the bundle axis
    tip_particles_flat = []
    for (top_end, bottom_start) in helix_tip_idx:
        tip_particles_flat.extend([top_end, bottom_start])
    tip_centroid = np.mean([all_positions[i] for i in tip_particles_flat], axis=0)
    root_particles_flat = []
    for (top_start, bottom_end) in helix_root_idx:
        root_particles_flat.extend([top_start, bottom_end])
    root_centroid = np.mean([all_positions[i] for i in root_particles_flat], axis=0)
    axis = tip_centroid - root_centroid
    axis /= np.linalg.norm(axis)
    # pick a lateral direction perpendicular to the bundle axis (x-ish)
    lateral = np.array([1.0, 0.0, 0.0])
    lateral -= axis * np.dot(axis, lateral)
    lateral /= np.linalg.norm(lateral)

    bias_target = tip_centroid + lateral * (bias_nm / OXDNA_LENGTH_NM)

    with open("ext.dat", "w") as out:
        for (top_start, bottom_end) in helix_root_idx:
            for idx in (top_start, bottom_end):
                x, y, z = all_positions[idx]
                out.write(f"{{\ntype = trap\nparticle = {idx}\nstiff = {trap_stiff}\n"
                          f"pos0 = {x}, {y}, {z}\nrate = 0.\ndir = 0.,0.,0.\n}}\n\n")

        def crosslink_ring(ring_idx_list):
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

        crosslink_ring(helix_tip_idx)

        # bias trap: translates each tip particle's OWN natural position
        # by the bias vector (not toward one shared point -- that would
        # fight the ring crosslinks, which enforce nonzero hexagon
        # spacing, and collapsing 12 points to one shared target breaks
        # bonds). Rigidly shifting the whole ring by the bias vector
        # preserves ring geometry while moving its mean position --
        # proxy for a strand-displacement-driven hinge lock.
        bias_vec = np.array(bias_target) - np.array(tip_centroid)
        lock_set = tip_particles_flat if n_lock_particles == 0 else tip_particles_flat[:n_lock_particles]
        for idx in lock_set:
            natural = np.array(all_positions[idx])
            target = natural + bias_vec
            x, y, z = target
            out.write(f"{{\ntype = trap\nparticle = {idx}\nstiff = {trap_stiff}\n"
                      f"pos0 = {x}, {y}, {z}\nrate = 0.\ndir = 0.,0.,0.\n}}\n\n")
        print(f"Lock applied to {len(lock_set)}/{len(tip_particles_flat)} tip particles: {lock_set}")

    with open("tip_particles.txt", "w") as out:
        for h in range(n_helix):
            out.write(f"{helix_tip_idx[h][0]} {helix_tip_idx[h][1]}\n")

    print(f"Built {n_helix}-helix hinge-lever, bp={bp}, bias={bias_nm}nm")
    print(f"Natural tip centroid: {tip_centroid}")
    print(f"Biased target: {bias_target}")
    print(f"Box side: {box_side:.2f} units")


if __name__ == "__main__":
    main()
