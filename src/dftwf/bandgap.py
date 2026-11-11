from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

_HOMO_LUMO_RE = re.compile(
    r"highest occupied, lowest unoccupied level \(ev\):\s*"
    r"(-?\d+\.\d+)\s+(-?\d+\.\d+)",
    re.I,
)
_HOMO_ONLY_RE = re.compile(r"highest occupied level \(ev\):\s*(-?\d+\.\d+)", re.I)
_NELEC_RE = re.compile(r"number of electrons\s*=\s*(\d+(?:\.\d+)?)", re.I)


@dataclass
class BandGap:
    """Result of a band gap analysis (energies in eV)."""

    vbm_ev: float
    cbm_ev: float
    gap_ev: float
    direct_gap_ev: float | None = None   # smallest Ec(k) - Ev(k), path only
    is_direct: bool | None = None        # None when not determined
    k_vbm: float | None = None           # path coordinate of the VBM
    k_cbm: float | None = None           # path coordinate of the CBM
    is_metal: bool = False


def get_number_of_electrons(path):
    """Return the number of electrons printed in a pw.x output."""
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    m = _NELEC_RE.search(text)
    if m is None:
        raise ValueError(f"'number of electrons' not found in {path}")
    return float(m.group(1))


def get_homo_lumo(path):
    """Return (VBM, CBM) in eV from the 'highest occupied, lowest unoccupied level' line."""
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    matches = _HOMO_LUMO_RE.findall(text)
    if matches:
        vbm, cbm = matches[-1]
        return float(vbm), float(cbm)
    if _HOMO_ONLY_RE.search(text):
        raise ValueError(
            f"Only the highest occupied level is printed in {path}; "
            "add empty bands (nbnd > nelec/2) to get the gap."
        )
    raise ValueError(
        f"No 'highest occupied, lowest unoccupied level' line in {path} "
        "(metal with smearing, or no empty bands in the run)."
    )


def gap_from_output(path):
    """Gap from the k-points of a pw.x run (scf/nscf with empty bands)."""
    vbm, cbm = get_homo_lumo(path)
    gap = cbm - vbm
    return BandGap(vbm_ev=vbm, cbm_ev=cbm, gap_ev=max(gap, 0.0), is_metal=gap <= 0)


def read_band_energies(path):
    """Read a bands.x gnuplot file. Returns (k, energies) with energies[nbands, nk]."""
    bands, current = [], []
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    for lineno, line in enumerate(text.splitlines(), start=1):
        parts = line.split()
        if not parts:
            if current:
                bands.append(current)
                current = []
            continue
        if len(parts) != 2:
            raise ValueError(
                f"{path}, line {lineno}: expected 2 columns (k, E), got {len(parts)}"
            )
        current.append((float(parts[0]), float(parts[1])))
    if current:
        bands.append(current)
    if not bands:
        raise ValueError(f"No band data found in {path}")
    if len({len(b) for b in bands}) != 1:
        raise ValueError(f"Bands in {path} have different numbers of k-points")
    arr = np.array(bands, dtype=float)  # (nbands, nk, 2)
    return arr[0, :, 0], arr[:, :, 1]


def _n_occupied(nelec):
    n = round(nelec)
    if abs(nelec - n) > 1e-6 or n % 2 != 0:
        raise ValueError(
            f"nelec = {nelec:g} is odd or non-integer: the number of occupied bands "
            "is ambiguous (metal or spin-polarized system?)."
        )
    return n // 2


def gap_from_bands(bands_file, nelec, direct_tol=1e-3):
    """Gap along the k-path of a bands.dat.gnu file (non-spin-polarized).

    Occupied bands = nelec / 2. Only valid along the path used: band extrema
    that lie off the path are missed.
    """
    n_occ = _n_occupied(nelec)
    k, e = read_band_energies(bands_file)
    if e.shape[0] <= n_occ:
        raise ValueError(
            f"Need more than {n_occ} bands to get a gap, found {e.shape[0]} in {bands_file}"
        )
    ev, ec = e[n_occ - 1], e[n_occ]
    iv, ic = int(np.argmax(ev)), int(np.argmin(ec))
    vbm, cbm = float(ev[iv]), float(ec[ic])
    gap = cbm - vbm
    direct_gap = float(np.min(ec - ev))
    metal = gap <= 0
    return BandGap(
        vbm_ev=vbm,
        cbm_ev=cbm,
        gap_ev=max(gap, 0.0),
        direct_gap_ev=direct_gap,
        is_direct=None if metal else (direct_gap - gap) < direct_tol,
        k_vbm=float(k[iv]),
        k_cbm=float(k[ic]),
        is_metal=metal,
    )


def analyze_gap(output, bands_file=None):
    """Return (gap_from_output_or_None, gap_from_bands_or_None)."""
    try:
        from_output = gap_from_output(output)
    except ValueError:
        from_output = None
    from_bands = None
    if bands_file:
        from_bands = gap_from_bands(bands_file, get_number_of_electrons(output))
    if from_output is None and from_bands is None:
        raise ValueError(
            f"No HOMO/LUMO line in {output} and no --bands file given; nothing to analyse."
        )
    return from_output, from_bands


def format_gap(label, g):
    """Human-readable multi-line summary of a BandGap."""
    lines = [label, f"  VBM = {g.vbm_ev:.4f} eV, CBM = {g.cbm_ev:.4f} eV"]
    if g.is_metal:
        lines.append("  Gap = 0 (bands overlap: metallic)")
    else:
        kind = ""
        if g.is_direct is not None:
            kind = " (direct)" if g.is_direct else " (indirect)"
        lines.append(f"  Gap = {g.gap_ev:.4f} eV{kind}")
    if g.direct_gap_ev is not None:
        lines.append(f"  Smallest direct gap = {g.direct_gap_ev:.4f} eV")
    if g.k_vbm is not None and g.k_cbm is not None:
        lines.append(
            f"  k(VBM) = {g.k_vbm:.4f}, k(CBM) = {g.k_cbm:.4f} (path coordinate)"
        )
    return "\n".join(lines)
