# oxDNA Bundle Toolkit

Tools for building and analyzing multi-helix DNA duplex bundles in
[oxDNA](https://github.com/lorenzo-rovigatti/oxDNA), a coarse-grained
molecular dynamics engine for DNA/RNA nanostructures. Built while
researching candidate mechanical elements (rigid cantilevers and
bistable hinge-levers) for a nanotechnology design project targeting
the [Foresight Institute Feynman Grand Prize](https://foresight.org/feynman-grand-prize/).

## What's here

- **`build_bundle.py`** — builds an N-helix DNA duplex bundle: a ring
  of parallel double helices at realistic DNA-nanotech spacing (2.5nm
  helix-to-helix), root ends individually anchored, tip ends
  crosslinked into a rigid ring via oxDNA external forces (an
  idealized proxy for real covalent crossovers — see the module
  docstring for why, and its limitations). Parametrized by base-pair
  length, helix count, number of interior crosslink rings, and trap
  stiffness — used to study how bundle geometry affects mechanical
  rigidity.
- **`build_hinge_lever.py`** — extends the bundle builder with a
  controllable bias force that displaces the tip toward an offset
  target, modeling a strand-displacement-triggered mechanical lock
  (e.g. a TMSD-driven DNA hinge settling into one of two states).
  Supports a configurable number of "lock particles" to compare an
  idealized full lock against a more realistic single-attachment-point
  lock.
- **`analyze_bundle_trajectory.py`** / **`analyze_cantilever_trajectory.py`**
  — parse an oxDNA trajectory and compute the tip's RMS positional
  fluctuation and the implied effective transverse stiffness (via the
  equipartition theorem, k = k_B·T / ⟨x²⟩).
- **`constants.py`** — shared physical constants and DNA-nanotech
  material-property reference values used across the scripts.

## Requirements

- Python 3.9+, NumPy
- [oxDNA](https://github.com/lorenzo-rovigatti/oxDNA) built and on
  your `PATH` (or reference its binary directly) to actually run the
  generated structures — these scripts produce oxDNA input files
  (`generated.top`, `generated.dat`, `ext.dat`) and analyze its output
  trajectories; they don't run the simulation themselves.

## Usage

```bash
# Build a 6-helix, 10bp bundle (default trap stiffness)
python3 build_bundle.py 10 6

# Build a 12-helix, 20bp bundle with 1 extra interior crosslink ring
python3 build_bundle.py 20 12 1

# Build a hinge-lever: 6-helix, 10bp arm, biased 3nm laterally,
# realistic 2-particle (single attachment point) lock
python3 build_hinge_lever.py 10 6 3.0 5.0 2

# Then run oxDNA on the generated files (see oxDNA's own docs for
# input-file / MD-parameter setup), and analyze the resulting
# trajectory:
python3 analyze_bundle_trajectory.py trajectory.dat tip_particles.txt
```

## Background

These tools were built to test a specific question: does DNA-origami
bundling close the gap between a floppy single duplex (thermally
limited to several nm of positional noise) and the sub-nanometer
rigidity needed for precision nanomechanical elements? Short answer:
bundling and shortening both help substantially (single duplex →
~3,100x too soft → 6-42 helix bundles bring that down to single-digit
multiples), but a real, quantified diagnostic (comparing idealized vs.
realistic lock/crosslink stiffness) shows a meaningful fraction of the
remaining softness lives in the anchor/joint points, not the beam
material itself — a finding independently consistent with known DNA
nanotechnology literature on crossover/joint flexibility.

## License

MIT — see `LICENSE`.
