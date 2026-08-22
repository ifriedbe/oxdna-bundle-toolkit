"""Physical constants and real material property estimates used across
the simulation scripts. Values are SI unless noted, sourced from
standard physical constants (CODATA) and commonly cited nanomechanics
literature ranges. Where a value is a rough order-of-magnitude estimate
rather than a precise measurement, it's noted as such.
"""

import numpy as np

kB = 1.380649e-23          # Boltzmann constant, J/K
amu = 1.66053906660e-27    # atomic mass unit, kg
angstrom = 1e-10           # m
nm = 1e-9                  # m

T_ROOM = 300.0              # K
T_LIQUID_N2 = 77.0          # K
T_LIQUID_HE = 4.2           # K

# --- Candidate actuator stiffness estimates (N/m) ---
# These are order-of-magnitude figures pulled from published nanomechanics
# work, used only to bound plausible design space, not exact predictions
# for any specific device.
STIFFNESS_DNA_ORIGAMI = 1e-3      # ~pN/nm scale hinges/springs (soft, entropic)
STIFFNESS_CARBON_NANOTUBE = 1.0    # order N/m for a short CNT cantilever segment
STIFFNESS_DIAMONDOID_SMALL = 50.0  # small rigid diamondoid strut, order N/m
STIFFNESS_DIAMOND_BULK_EQUIV = 500.0  # very short, very stiff covalent strut

# --- Candidate arm/tip effective mass estimates (kg) ---
# A "few thousand to a few hundred thousand atoms" manipulator tip/arm.
MASS_SMALL_CLUSTER_1K_ATOMS = 1000 * 12 * amu     # ~1000 carbon atoms
MASS_MED_CLUSTER_100K_ATOMS = 100_000 * 12 * amu  # ~100,000 carbon atoms
MASS_LARGE_CLUSTER_1M_ATOMS = 1_000_000 * 12 * amu

# --- Lennard-Jones parameters (illustrative, carbon-carbon-like) ---
# epsilon in Joules, sigma in meters. These are generic sp3-carbon-scale
# LJ parameters used for a first-principles vibrational-frequency sanity
# check, not a validated force field for any specific candidate material.
LJ_EPSILON_CC = 0.00284 * 1.602176634e-19  # ~0.00284 eV -> J (graphite-like vdW well depth order)
LJ_SIGMA_CC = 3.4 * angstrom                # ~3.4 Angstrom (graphite interlayer-like)

# Covalent C-C bond stiffness (from diamond/graphene literature), used as
# an upper bound for how stiff a real atomic-scale strut can be.
COVALENT_CC_BOND_STIFFNESS = 500.0  # N/m, commonly cited order for a single C-C bond stretch mode
