"""Analyze the 6-helix bundle trajectory: RMS fluctuation of the tip
cross-section's CENTROID (the position that would encode the output
bit, analogous to the single-duplex tip in analyze_cantilever_
trajectory.py), and the implied effective stiffness for direct
comparison against the single-duplex result (1.3e-4 N/m, 3100x too
soft vs. the 0.41 N/m target).
"""

import sys
import numpy as np

OXDNA_LENGTH_NM = 0.8518


def parse_trajectory(path, particle_indices):
    traj = {i: [] for i in particle_indices}
    with open(path) as f:
        lines = f.readlines()
    i, n = 0, len(lines)
    while i < n:
        if lines[i].startswith("t ="):
            i += 3
            particle_idx = 0
            while i < n and not lines[i].startswith("t ="):
                if particle_idx in traj:
                    parts = lines[i].split()
                    traj[particle_idx].append([float(parts[0]), float(parts[1]), float(parts[2])])
                particle_idx += 1
                i += 1
        else:
            i += 1
    return {k: np.array(v) for k, v in traj.items()}


def main():
    traj_path = sys.argv[1] if len(sys.argv) > 1 else "trajectory.dat"
    tip_file = sys.argv[2] if len(sys.argv) > 2 else "tip_particles.txt"

    tip_particle_ids = []
    with open(tip_file) as f:
        for line in f:
            a, b = line.split()
            tip_particle_ids.extend([int(a), int(b)])

    traj = parse_trajectory(traj_path, tip_particle_ids)
    n_frames = len(traj[tip_particle_ids[0]])
    print(f"Frames captured: {n_frames}")
    print(f"Tip particles (12, 2 per helix x 6 helices): {tip_particle_ids}")

    # centroid of all 12 tip particles per frame = the bundle's tip position
    stacked = np.stack([traj[i] for i in tip_particle_ids], axis=0)  # (12, frames, 3)
    centroid = stacked.mean(axis=0)  # (frames, 3)

    mean = centroid.mean(axis=0)
    disp = centroid - mean
    rms_nm = np.sqrt((disp ** 2).sum(axis=1).mean()) * OXDNA_LENGTH_NM
    max_nm = np.sqrt((disp ** 2).sum(axis=1)).max() * OXDNA_LENGTH_NM

    print(f"\nbundle tip centroid: RMS fluctuation = {rms_nm:.4f} nm, max excursion = {max_nm:.4f} nm")

    kB = 1.380649e-23
    T = 300.0
    x2_m2 = (rms_nm * 1e-9) ** 2
    k_eff = kB * T / x2_m2
    print(f"implied effective transverse stiffness: {k_eff:.3e} N/m")
    print(f"vs. 0.41 N/m target: {'MEETS/EXCEEDS' if k_eff >= 0.41 else 'still short by'} "
          f"{'' if k_eff >= 0.41 else f'{0.41 / k_eff:.0f}x'}")

    single_duplex_k = 1.322e-4
    print(f"\nvs. single-duplex result (1.322e-4 N/m): "
          f"{k_eff / single_duplex_k:.1f}x stiffer")

    print(f"\n[Reminder: this bundle uses idealized rigid trap crosslinks at the")
    print(f" tip, not real covalent crossovers -- an upper-bound proxy. See")
    print(f" build_bundle.py docstring.]")


if __name__ == "__main__":
    main()
