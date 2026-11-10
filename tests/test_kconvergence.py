import shutil
from pathlib import Path

import matplotlib
import pytest

matplotlib.use("Agg")

from dftwf.convergence import (
    analyze_kpoints,
    collect_kpoint_results,
    find_kconverged,
    get_kgrid_from_input,
    get_nk_irreducible,
)

ROOT = Path(__file__).resolve().parent.parent
KCONV = ROOT / "examples" / "si" / "kconvergence"


def test_parse_kgrid_and_nk_ibz():
    assert get_kgrid_from_input(str(KCONV / "scf_k12.in")) == (12, 12, 12)
    assert get_kgrid_from_input(str(KCONV / "scf_k2.in")) == (2, 2, 2)
    assert get_nk_irreducible(str(KCONV / "scf_k12.out")) == 72
    assert get_nk_irreducible(str(KCONV / "scf_k2.out")) == 3


def test_kgrid_missing_block_raises(tmp_path):
    inp = tmp_path / "bad.in"
    inp.write_text("&CONTROL\n calculation = 'scf'\n/\n")
    with pytest.raises(ValueError):
        get_kgrid_from_input(str(inp))


def test_collect_sorted_valid_and_deltas():
    points = collect_kpoint_results(KCONV)
    assert [p.grid for p in points] == [
        (2, 2, 2), (4, 4, 4), (6, 6, 6), (8, 8, 8), (12, 12, 12),
    ]
    assert [p.nk_irreducible for p in points] == [3, 8, 16, 29, 72]
    assert all(p.valid for p in points)
    assert points[0].delta_mev_per_atom is None
    expected = [812.7785, 15.4438, 0.0896, 0.0002]
    for point, value in zip(points[1:], expected):
        assert point.delta_mev_per_atom == pytest.approx(value, abs=1e-3)


@pytest.mark.parametrize(
    "threshold, expected",
    [(20.0, (6, 6, 6)), (1.0, (8, 8, 8)), (0.1, (8, 8, 8)), (1e-5, None)],
)
def test_find_kconverged(threshold, expected):
    points = collect_kpoint_results(KCONV)
    result = find_kconverged(points, threshold_mev=threshold)
    assert (None if result is None else result.grid) == expected


def test_unfinished_job_is_invalid(tmp_path):
    for src in KCONV.glob("scf_k*.in"):
        shutil.copy(src, tmp_path / src.name)
    for src in KCONV.glob("scf_k*.out"):
        text = src.read_text()
        if src.name == "scf_k12.out":
            text = "\n".join(
                line for line in text.splitlines() if "JOB DONE" not in line
            )
        (tmp_path / src.name).write_text(text)

    points = collect_kpoint_results(tmp_path)
    invalid = [p for p in points if not p.valid]
    assert [p.file for p in invalid] == ["scf_k12.out"]
    assert find_kconverged(points, threshold_mev=0.1).grid == (8, 8, 8)


def test_missing_input_file_is_invalid(tmp_path):
    for src in KCONV.glob("scf_k*.out"):
        shutil.copy(src, tmp_path / src.name)
    for src in KCONV.glob("scf_k*.in"):
        if src.name != "scf_k12.in":
            shutil.copy(src, tmp_path / src.name)

    points = collect_kpoint_results(tmp_path)
    invalid = [p for p in points if not p.valid]
    assert [p.file for p in invalid] == ["scf_k12.out"]
    assert invalid[0].grid is None


def test_empty_folder_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        collect_kpoint_results(tmp_path)


def test_analyze_writes_csv_and_plot(tmp_path):
    csv_path = tmp_path / "kconv.csv"
    png_path = tmp_path / "kconv.png"
    points, converged = analyze_kpoints(
        KCONV, threshold_mev=0.1, csv_path=csv_path, plot_path=png_path
    )
    assert converged.grid == (8, 8, 8)
    assert png_path.exists() and png_path.stat().st_size > 0
    assert len(csv_path.read_text().strip().splitlines()) == 6  # header + 5 titik

def test_cli_kconverge(capsys):
    from dftwf.cli import main

    code = main(["kconverge", str(KCONV), "--threshold", "0.1"])
    out = capsys.readouterr().out
    assert code == 0
    assert "12x12x12" in out
    assert "Converged at k grid = 8x8x8" in out

def test_cli_kconverge_empty_folder(tmp_path, capsys):
    from dftwf.cli import main

    assert main(["kconverge", str(tmp_path)]) == 1
    assert "Error" in capsys.readouterr().err
