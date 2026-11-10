from __future__ import annotations

import csv
import re
import warnings
from dataclasses import dataclass
from pathlib import Path

from .parsers import (
    get_ecutwfc,
    get_number_of_atoms,
    get_total_energy,
    is_converged,
    is_job_done,
)

RY_TO_MEV = 13605.693122994  # 1 Ry in meV


@dataclass
class ConvergencePoint:
    """One pw.x run in a convergence series."""

    file: str
    ecutwfc: float | None
    n_atoms: int | None
    energy_ry: float | None
    valid: bool
    delta_mev_per_atom: float | None = None


def _safe(func, path):
    """Call a parser and return None if it raises."""
    try:
        return func(path)
    except Exception:
        return None


def _compute_deltas(valid_points):
    """Fill delta_mev_per_atom: |E_i - E_(i-1)| in meV/atom (points sorted by ecut)."""
    for prev, cur in zip(valid_points, valid_points[1:]):
        cur.delta_mev_per_atom = (
            abs(cur.energy_ry - prev.energy_ry) * RY_TO_MEV / cur.n_atoms
        )


def collect_ecut_results(folder, pattern="*.out"):
    """Read all pw.x outputs in a folder and return a list of ConvergencePoint.

    Valid points (job finished, SCF converged, all values found) come first,
    sorted by ecutwfc. Invalid points follow and are excluded from the analysis.
    """
    paths = sorted(Path(folder).glob(pattern))
    if not paths:
        raise FileNotFoundError(f"No files matching '{pattern}' in {folder}")

    points = []
    for path in paths:
        p = str(path)
        ecut = _safe(get_ecutwfc, p)
        n_atoms = _safe(get_number_of_atoms, p)
        energy = _safe(get_total_energy, p)
        ok = (
            bool(_safe(is_job_done, p))
            and bool(_safe(is_converged, p))
            and ecut is not None
            and n_atoms is not None
            and energy is not None
        )
        points.append(ConvergencePoint(path.name, ecut, n_atoms, energy, ok))

    valid = sorted((pt for pt in points if pt.valid), key=lambda pt: pt.ecutwfc)
    invalid = [pt for pt in points if not pt.valid]
    _compute_deltas(valid)
    return valid + invalid


def find_converged(points, threshold_mev=1.0):
    """Return the lowest ecutwfc after which every step is below the threshold.

    The step is |E(i) - E(i-1)| in meV/atom. Returns None if never reached.
    """
    valid = [pt for pt in points if pt.valid]
    deltas = [pt.delta_mev_per_atom for pt in valid]
    for i in range(1, len(valid)):
        if all(d < threshold_mev for d in deltas[i:]):
            return valid[i].ecutwfc
    return None


def write_csv(points, path):
    """Write the convergence table to a CSV file."""
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(
            ["file", "ecutwfc_ry", "n_atoms", "total_energy_ry",
             "delta_mev_per_atom", "valid"]
        )
        for pt in points:
            writer.writerow([
                pt.file,
                "" if pt.ecutwfc is None else pt.ecutwfc,
                "" if pt.n_atoms is None else pt.n_atoms,
                "" if pt.energy_ry is None else pt.energy_ry,
                "" if pt.delta_mev_per_atom is None else f"{pt.delta_mev_per_atom:.6f}",
                pt.valid,
            ])


def plot_convergence(points, output, threshold_mev=1.0, converged=None, title=None):
    """Plot E - E(last point) and the step size |dE| against ecutwfc."""
    import matplotlib.pyplot as plt

    valid = [pt for pt in points if pt.valid]
    if len(valid) < 2:
        raise ValueError("Need at least two valid points to plot convergence.")

    x = [pt.ecutwfc for pt in valid]
    ref = valid[-1]
    rel = [(pt.energy_ry - ref.energy_ry) * RY_TO_MEV / pt.n_atoms for pt in valid]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)

    ax1.plot(x, rel, "o-")
    ax1.axhline(0.0, color="gray", linestyle=":", linewidth=1)
    ax1.set_xlabel("ecutwfc (Ry)")
    ax1.set_ylabel("E − E(last point) (meV/atom)")

    ax2.semilogy(x[1:], [pt.delta_mev_per_atom for pt in valid[1:]], "s-")
    ax2.axhline(
        threshold_mev, color="gray", linestyle="--",
        label=f"threshold = {threshold_mev:g} meV/atom",
    )
    if converged is not None:
        ax2.axvline(
            converged, color="tab:green", linestyle=":",
            label=f"converged at {converged:g} Ry",
        )
    ax2.set_xlabel("ecutwfc (Ry)")
    ax2.set_ylabel("|ΔE| between steps (meV/atom)")
    ax2.legend()

    if title:
        fig.suptitle(title)
    fig.savefig(output, dpi=300)
    plt.close(fig)


def analyze_ecut(folder, threshold_mev=1.0, csv_path=None, plot_path=None,
                 pattern="*.out", title=None):
    """Collect results, find the converged cutoff, optionally write CSV and plot.

    Returns (points, converged_ecutwfc_or_None).
    """
    points = collect_ecut_results(folder, pattern=pattern)
    converged = find_converged(points, threshold_mev=threshold_mev)
    if csv_path:
        write_csv(points, csv_path)
    if plot_path:
        plot_convergence(points, plot_path, threshold_mev=threshold_mev,
                         converged=converged, title=title)
    return points, converged

_KAUTO_RE = re.compile(r"^\s*K_POINTS\s*[{(]?\s*automatic\s*[})]?\s*(!.*)?$", re.I)
_NK_IBZ_RE = re.compile(r"number of k points=\s*(\d+)")


@dataclass
class KpointPoint:
    """One pw.x run in a k-point convergence series."""

    file: str
    grid: tuple[int, int, int] | None
    nk_grid: int | None          # n1*n2*n3 (sorting key)
    nk_irreducible: int | None   # number of k-points in the IBZ (from output)
    ecutwfc: float | None
    n_atoms: int | None
    energy_ry: float | None
    valid: bool
    delta_mev_per_atom: float | None = None

    @property
    def grid_label(self):
        return "-" if self.grid is None else "x".join(map(str, self.grid))


def get_kgrid_from_input(path):
    """Return (n1, n2, n3) from a 'K_POINTS automatic' block in a pw.x input."""
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    for i, line in enumerate(lines):
        if _KAUTO_RE.match(line):
            for nxt in lines[i + 1:]:
                if nxt.strip():
                    parts = nxt.split()
                    if len(parts) < 3:
                        break
                    return int(parts[0]), int(parts[1]), int(parts[2])
            break
    raise ValueError(f"No 'K_POINTS automatic' grid found in {path}")


def get_nk_irreducible(path):
    """Return the number of k-points in the IBZ as printed in a pw.x output."""
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    m = _NK_IBZ_RE.search(text)
    if m is None:
        raise ValueError(f"'number of k points' not found in {path}")
    return int(m.group(1))


def collect_kpoint_results(folder, pattern="*.out", input_suffix=".in"):
    """Read pw.x outputs (+ matching inputs) and return a list of KpointPoint.

    The k grid is read from the input file with the same stem as the output.
    Valid points come first, sorted by total grid points; invalid points follow.
    """
    paths = sorted(Path(folder).glob(pattern))
    if not paths:
        raise FileNotFoundError(f"No files matching '{pattern}' in {folder}")

    points = []
    for path in paths:
        p = str(path)
        inp = path.with_suffix(input_suffix)
        grid = _safe(get_kgrid_from_input, str(inp)) if inp.exists() else None
        nk_grid = None if grid is None else grid[0] * grid[1] * grid[2]
        nk_ibz = _safe(get_nk_irreducible, p)
        ecut = _safe(get_ecutwfc, p)
        n_atoms = _safe(get_number_of_atoms, p)
        energy = _safe(get_total_energy, p)
        ok = (
            bool(_safe(is_job_done, p))
            and bool(_safe(is_converged, p))
            and grid is not None
            and n_atoms is not None
            and energy is not None
        )
        points.append(
            KpointPoint(path.name, grid, nk_grid, nk_ibz, ecut, n_atoms, energy, ok)
        )

    valid = sorted((pt for pt in points if pt.valid),
                   key=lambda pt: (pt.nk_grid, pt.grid))
    invalid = [pt for pt in points if not pt.valid]

    ecuts = {pt.ecutwfc for pt in valid if pt.ecutwfc is not None}
    if len(ecuts) > 1:
        warnings.warn(
            f"ecutwfc differs between runs ({sorted(ecuts)} Ry); "
            "a k-point test should keep ecutwfc fixed."
        )

    _compute_deltas(valid)
    return valid + invalid


def find_kconverged(points, threshold_mev=1.0):
    """Return the first KpointPoint after which every step is below the threshold."""
    valid = [pt for pt in points if pt.valid]
    deltas = [pt.delta_mev_per_atom for pt in valid]
    for i in range(1, len(valid)):
        if all(d < threshold_mev for d in deltas[i:]):
            return valid[i]
    return None


def write_kpoint_csv(points, path):
    """Write the k-point convergence table to a CSV file."""
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(
            ["file", "kgrid", "nk_grid", "nk_irreducible", "ecutwfc_ry",
             "n_atoms", "total_energy_ry", "delta_mev_per_atom", "valid"]
        )
        for pt in points:
            writer.writerow([
                pt.file,
                pt.grid_label,
                "" if pt.nk_grid is None else pt.nk_grid,
                "" if pt.nk_irreducible is None else pt.nk_irreducible,
                "" if pt.ecutwfc is None else pt.ecutwfc,
                "" if pt.n_atoms is None else pt.n_atoms,
                "" if pt.energy_ry is None else pt.energy_ry,
                "" if pt.delta_mev_per_atom is None else f"{pt.delta_mev_per_atom:.6f}",
                pt.valid,
            ])


def plot_kconvergence(points, output, threshold_mev=1.0, converged=None, title=None):
    """Plot E - E(last point) and the step size |dE| against the k grid."""
    import matplotlib.pyplot as plt

    valid = [pt for pt in points if pt.valid]
    if len(valid) < 2:
        raise ValueError("Need at least two valid points to plot convergence.")

    idx = list(range(len(valid)))
    labels = [pt.grid_label for pt in valid]
    ref = valid[-1]
    rel = [(pt.energy_ry - ref.energy_ry) * RY_TO_MEV / pt.n_atoms for pt in valid]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)

    ax1.plot(idx, rel, "o-")
    ax1.axhline(0.0, color="gray", linestyle=":", linewidth=1)
    ax1.set_xticks(idx, labels, rotation=45)
    ax1.set_xlabel("k-point grid")
    ax1.set_ylabel("E − E(last point) (meV/atom)")

    ax2.semilogy(idx[1:], [pt.delta_mev_per_atom for pt in valid[1:]], "s-")
    ax2.axhline(
        threshold_mev, color="gray", linestyle="--",
        label=f"threshold = {threshold_mev:g} meV/atom",
    )
    if converged is not None:
        ax2.axvline(
            valid.index(converged), color="tab:green", linestyle=":",
            label=f"converged at {converged.grid_label}",
        )
    ax2.set_xticks(idx, labels, rotation=45)
    ax2.set_xlabel("k-point grid")
    ax2.set_ylabel("|ΔE| between steps (meV/atom)")
    ax2.legend()

    if title:
        fig.suptitle(title)
    fig.savefig(output, dpi=300)
    plt.close(fig)


def analyze_kpoints(folder, threshold_mev=1.0, csv_path=None, plot_path=None,
                    pattern="*.out", title=None):
    """Collect results, find the converged k grid, optionally write CSV and plot.

    Returns (points, converged_KpointPoint_or_None).
    """
    points = collect_kpoint_results(folder, pattern=pattern)
    converged = find_kconverged(points, threshold_mev=threshold_mev)
    if csv_path:
        write_kpoint_csv(points, csv_path)
    if plot_path:
        plot_kconvergence(points, plot_path, threshold_mev=threshold_mev,
                          converged=converged, title=title)
    return points, converged
