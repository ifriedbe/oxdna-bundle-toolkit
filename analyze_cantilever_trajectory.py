"""Parse an oxDNA trajectory.dat and compute the RMS positional
fluctuation of the free (distal) end of an anchored DNA-duplex
'cantilever', as a real-MD cross-check against an analytical
thermal-noise estimate (k = k_B*T / <x^2>, e.g. k ~= 0.41 N/m needed
for 0.1nm positional accuracy at 300K in this project's own use case).

oxDNA length unit = 0.8518 nm (standard oxDNA unit conversion).
"""

import sys
import numpy as np

OXDNA_LENGTH_NM = 0.8518

def parse_trajectory(path, particle_indices):
    """Return dict: particle_index -> Nx3 array of positions (oxDNA units)."""
    traj = {i: [] for i in particle_indices}
    with open(path) as f:
        lines = f.readlines()
    i = 0
    n = len(lines)
    while i < n:
        if lines[i].startswith("t ="):
            # header: t=, b=, E=
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
    path = sys.argv[1] if len(sys.argv) > 1 else "trajectory.dat"
    tip_particles = (19, 20)
    root_particles = (0, 39)

    traj = parse_trajectory(path, tip_particles + root_particles)

    print(f"Frames captured: {len(traj[19])}")

    for label, idx in (("tip (19)", 19), ("tip-pair (20)", 20),
                        ("root (0, should be ~fixed)", 0), ("root-pair (39, should be ~fixed)", 39)):
        pos = traj[idx]
        mean = pos.mean(axis=0)
        disp = pos - mean
        rms_nm = np.sqrt((disp ** 2).sum(axis=1).mean()) * OXDNA_LENGTH_NM
        max_nm = np.sqrt((disp ** 2).sum(axis=1)).max() * OXDNA_LENGTH_NM
        print(f"{label}: RMS fluctuation = {rms_nm:.4f} nm, max excursion = {max_nm:.4f} nm")

    # average tip position (midpoint of the two paired tip bases) as the
    # effective 'cantilever tip' whose height would encode the output bit
    tip_mid = (traj[19] + traj[20]) / 2.0
    mean = tip_mid.mean(axis=0)
    disp = tip_mid - mean
    rms_nm = np.sqrt((disp ** 2).sum(axis=1).mean()) * OXDNA_LENGTH_NM
    print(f"\ntip midpoint RMS fluctuation: {rms_nm:.4f} nm")

    # Effective transverse stiffness implied by equipartition, k = kB*T / <x^2>.
    # The 0.41 N/m reference value below is this project's own example target
    # (the stiffness needed for 0.1nm positional accuracy at 300K); swap in
    # whatever target applies to your own use case.
    kB = 1.380649e-23
    T = 300.0
    x2_m2 = (rms_nm * 1e-9) ** 2
    k_eff = kB * T / x2_m2
    print(f"implied effective transverse stiffness: {k_eff:.3e} N/m")
    print(f"(vs. example 0.41 N/m target -- {0.41 / k_eff:.0f}x too soft)")
    print(f"\n[note: 0.41 N/m corresponds to 0.1nm RMS positional accuracy at 300K;")
    print(f" adjust to whatever accuracy/temperature target applies to your design.]")


if __name__ == "__main__":
    main()
