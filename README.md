# qe-workflow-tools

[![tests](https://github.com/otnielxs/qe-workflow-tools/actions/workflows/tests.yml/badge.svg)](https://github.com/otnielxs/qe-workflow-tools/actions/workflows/tests.yml)

Python tools for Quantum ESPRESSO (QE) workflows: parse `pw.x`
outputs, read band structure and DOS data, make publication-ready plots, and convergency analysis.
Dependencies are kept minimal (`numpy`, `matplotlib`).


## Features

- Parse `pw.x` output: total energy (Ry), Fermi energy (eV), SCF convergence status, job completion marker
- Read band structure (`bands.dat.gnu`) and DOS data
- Read high-symmetry k-points and labels from `crystal_b` k-path inputs
- Combined band structure + DOS plot with `plot_band_dos`
- Cutoff-energy convergence analysis (CSV table, plot, `dftwf converge` command)
- Unit tests built on real QE outputs
- k-point convergence analysis (CSV table, plot, `dftwf kconverge` command)
- Band gap analysis (direct/indirect, from the pw.x output or `bands.dat.gnu`, `dftwf gap` command)

## Installation

````bash
git clone https://github.com/otnielxs/qe-workflow-tools.git
cd qe-workflow-tools
python -m pip install -e ".[dev]"
````

Requires Python 3.10 or newer.

## Quick start

Everything below uses the silicon example in `examples/si/`.

````python
from dftwf.parsers import get_total_energy, get_fermi_energy, is_converged, is_job_done
from dftwf.plotting import plot_band_dos

scf = "examples/si/scf.out"
print(get_total_energy(scf), "Ry")
print(get_fermi_energy(scf), "eV")
print(is_converged(scf), is_job_done(scf))

plot_band_dos(
    bands_file="examples/si/bands.dat.gnu",
    dos_file="examples/si/dos.dat",
    fermi_energy=get_fermi_energy(scf),
    kpath_output="examples/si/band.out",
    kpath_input="examples/si/bands.in",
    output="silicon_bands_dos.png",
    energy_range=(-12, 10),
    title="Electronic Properties of Silicon (Si)",
)
````

Or run `python examples/si/run_example.py`, or open `examples/si/demo.ipynb`.

> The pseudopotential for Si is not included. Download it from the
> Quantum ESPRESSO pseudopotential library to rerun the calculations
> from the `.in` files.

## Convergence test (cutoff energy)

Run `pw.x` for several `ecutwfc` values, put the outputs in one folder, then:

```bash
dftwf converge examples/si/convergence --threshold 0.1 \
    --csv convergence.csv --plot convergence.png
```
Output:

```text
file                  ecutwfc (Ry)            E (Ry)   dE (meV/atom)  status
scf_40.out                      40      -93.41380748               -  ok
scf_50.out                      50      -93.41390154          0.6399  ok
scf_60.out                      60      -93.41392143          0.1353  ok
scf_70.out                      70      -93.41393301          0.0788  ok
scf_80.out                      80      -93.41394229          0.0631  ok

Converged at ecutwfc = 70 Ry (threshold 0.1 meV/atom).
```

or from Python:

```python
from dftwf.convergence import analyze_ecut

points, converged = analyze_ecut(
    "examples/si/convergence",
    threshold_mev=0.1,
    csv_path="convergence.csv",
    plot_path="convergence.png",
)
print("Converged at ecutwfc =", converged, "Ry")
```

A cutoff counts as converged when every energy step from that point onward is
below the threshold (in meV/atom). Runs that did not finish or did not converge
are flagged as invalid and excluded from the analysis. The silicon data in `examples/si/` is for demonstration only; its smearing settings are not recommended production parameters.

## Convergence test (k-points)

Run `pw.x` for several k-point grids with everything else fixed (including
`ecutwfc` and, for metals, smearing and `degauss`). Put each output next to
its input with the same base name (`scf_k4.in` / `scf_k4.out`), then:

```bash
dftwf kconverge examples/si/kconvergence --threshold 0.1 \
    --csv kconvergence.csv --plot kconvergence.png
```

Output:

```
file                      k grid   nk(IBZ)            E (Ry)   dE (meV/atom)  status
scf_k2.out                 2x2x2         3      -93.29217499               -  ok
scf_k4.out                 4x4x4         8      -93.41165123        812.7785  ok
scf_k6.out                 6x6x6        16      -93.41392143         15.4438  ok
scf_k8.out                 8x8x8        29      -93.41393460          0.0896  ok
scf_k12.out             12x12x12        72      -93.41393457          0.0002  ok

Converged at k grid = 8x8x8 (threshold 0.1 meV/atom).
```

or from Python:

```python
from dftwf.convergence import analyze_kpoints

points, converged = analyze_kpoints(
    "examples/si/kconvergence",
    threshold_mev=0.1,
    csv_path="kconvergence.csv",
    plot_path="kconvergence.png",
)
print("Converged at k grid =", converged.grid_label)
```

The k grid is read from the `K_POINTS automatic` block of the input file
(pw.x outputs do not print it); the number of irreducible k-points is read
from the output. Runs are ordered by total grid points (n1·n2·n3). The
convergence criterion is the same as for the cutoff test. A warning is issued
if `ecutwfc` differs between runs. The silicon data in `examples/si/kconvergence/`
is for demonstration only: it uses `degauss = 0.2` Ry, which is far larger than
appropriate for a semiconductor and makes the energy converge unrealistically fast.

## Band gap

Two ways to get the gap from files you already have:

```bash
dftwf gap examples/si/nscf.out --bands examples/si/bands.dat.gnu
```

Output:

```
From pw.x output (k-points of that run):
  VBM = 6.2096 eV, CBM = 6.7817 eV
  Gap = 0.5721 eV

From band structure (valid along this k-path only):
  VBM = 6.1187 eV, CBM = 6.7317 eV
  Gap = 0.6130 eV (indirect)
  Smallest direct gap = 2.6016 eV
  k(VBM) = 0.0000, k(CBM) = 0.8500 (path coordinate)
```

- **From the pw.x output**: the `highest occupied, lowest unoccupied level`
  line that pw.x prints when the run includes empty bands. It only covers the
  k-points of that run.
- **From the band structure** (`bands.dat.gnu`): the number of occupied bands is
  `nelec/2`, with `nelec` read from the pw.x output. The gap is
  `min(CBM) - max(VBM)` along the k-path; it is called direct when it equals the
  smallest direct gap along the path. It is only valid along the path used:
  band extrema off the path are missed.

In the silicon example the two results come from different runs (nscf: fixed occupations, `nbnd = 20`, 24x24x24 grid; bands: smearing run), so their absolute VBM/CBM values are not directly comparable.

or from Python:

```python
from dftwf.bandgap import gap_from_bands, get_number_of_electrons

nelec = get_number_of_electrons("examples/si/nscf.out")
g = gap_from_bands("examples/si/bands.dat.gnu", nelec)
print(g.gap_ev, g.is_direct)
```

## Units

Total energies are returned in Rydberg (Ry) and Fermi energies in electronvolt (eV),
exactly as printed by QE.

## Project structure

````
src/dftwf/      parsers.py, kpath.py, plotting.py, convergence.py, bandgap.py, cli.py
tests/          unit tests and fixtures
examples/si/    silicon example (QE inputs/outputs, script, notebook)
````

## Tests

````bash
pytest -v
````
## Limitations

- Parsers are tested on Quantum ESPRESSO `pw.x` output (non-spin-polarized Si examples); other QE versions or output formats may need adjustments.
- Convergence analysis covers `ecutwfc` and k-point grids. Each k-point run needs its `.in` file next to the `.out`, and only `K_POINTS automatic` grids are supported (not explicit k-lists). All runs in a series should use the same shift, cutoff and smearing.
- Band gap analysis assumes a non-spin-polarized system with an even number of electrons. A metal is reported as zero gap when the bands overlap. The gap from `bands.dat.gnu` is only valid along the k-path used.
