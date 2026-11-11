from pathlib import Path

import pytest

from dftwf.bandgap import (
    analyze_gap,
    gap_from_bands,
    gap_from_output,
    get_homo_lumo,
    get_number_of_electrons,
    read_band_energies,
)
from dftwf.cli import main

ROOT = Path(__file__).resolve().parent.parent
SI = ROOT / "examples" / "si"
NSCF = SI / "nscf.out"
BANDS = SI / "bands.dat.gnu"


def _write_gnu(path, bands):
    lines = []
    for band in bands:
        for i, energy in enumerate(band):
            lines.append(f"{i * 0.1:.4f} {energy:.4f}")
        lines.append("")
    path.write_text("\n".join(lines) + "\n")


def test_number_of_electrons():
    assert get_number_of_electrons(str(NSCF)) == pytest.approx(8.0)


def test_homo_lumo_and_gap_from_output():
    vbm, cbm = get_homo_lumo(str(NSCF))
    assert vbm == pytest.approx(6.2096, abs=1e-4)
    assert cbm == pytest.approx(6.7817, abs=1e-4)
    g = gap_from_output(str(NSCF))
    assert g.gap_ev == pytest.approx(0.5721, abs=1e-4)
    assert not g.is_metal


def test_no_homo_lumo_line_raises(tmp_path):
    out = tmp_path / "x.out"
    out.write_text("     number of electrons       =         8.00\n")
    with pytest.raises(ValueError):
        get_homo_lumo(str(out))


def test_homo_only_line_raises(tmp_path):
    out = tmp_path / "x.out"
    out.write_text("     highest occupied level (ev):     6.2096\n")
    with pytest.raises(ValueError, match="empty bands"):
        get_homo_lumo(str(out))


def test_read_band_energies_shape():
    k, e = read_band_energies(str(BANDS))
    assert e.shape[1] == len(k)
    assert e.shape[0] > 4


def test_gap_from_bands_si_is_indirect():
    g = gap_from_bands(str(BANDS), 8.0)
    assert g.vbm_ev < g.cbm_ev
    assert g.gap_ev == pytest.approx(0.5721, abs=0.15)
    assert g.is_direct is False
    assert g.direct_gap_ev > g.gap_ev + 0.5
    assert not g.is_metal


def test_odd_electrons_raises():
    with pytest.raises(ValueError):
        gap_from_bands(str(BANDS), 7.0)


def test_too_few_bands_raises(tmp_path):
    f = tmp_path / "b.gnu"
    _write_gnu(f, [[0.0, 0.1], [1.0, 1.1]])
    with pytest.raises(ValueError, match="Need more than"):
        gap_from_bands(str(f), 8.0)


def test_bad_format_raises(tmp_path):
    f = tmp_path / "b.gnu"
    f.write_text("0.0 1.0 2.0\n")
    with pytest.raises(ValueError, match="2 columns"):
        read_band_energies(str(f))


def test_direct_gap_synthetic(tmp_path):
    f = tmp_path / "b.gnu"
    _write_gnu(f, [[-1.0, 0.0, -1.0], [1.0, 1.0, 1.0]])
    g = gap_from_bands(str(f), 2.0)
    assert g.gap_ev == pytest.approx(1.0)
    assert g.is_direct is True


def test_indirect_gap_synthetic(tmp_path):
    f = tmp_path / "b.gnu"
    _write_gnu(f, [[0.0, -1.0, -1.0], [2.0, 1.0, 2.0]])
    g = gap_from_bands(str(f), 2.0)
    assert g.gap_ev == pytest.approx(1.0)
    assert g.direct_gap_ev == pytest.approx(2.0)
    assert g.is_direct is False


def test_metal_synthetic(tmp_path):
    f = tmp_path / "b.gnu"
    _write_gnu(f, [[0.0, 1.0, 2.0], [0.5, 0.2, 3.0]])
    g = gap_from_bands(str(f), 2.0)
    assert g.is_metal
    assert g.gap_ev == 0.0
    assert g.is_direct is None


def test_analyze_gap_needs_something(tmp_path):
    out = tmp_path / "x.out"
    out.write_text("     number of electrons       =         8.00\n")
    with pytest.raises(ValueError):
        analyze_gap(str(out))


def test_cli_gap(capsys):
    code = main(["gap", str(NSCF), "--bands", str(BANDS)])
    out = capsys.readouterr().out
    assert code == 0
    assert "0.5721" in out
    assert "(indirect)" in out


def test_cli_gap_missing_file(tmp_path, capsys):
    assert main(["gap", str(tmp_path / "none.out")]) == 1
    assert "Error" in capsys.readouterr().err
