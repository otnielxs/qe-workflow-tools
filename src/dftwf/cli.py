from __future__ import annotations

import argparse
import sys

from .convergence import analyze_ecut, analyze_kpoints

def build_parser():
    parser = argparse.ArgumentParser(
        prog="dftwf",
        description="Lightweight tools for Quantum ESPRESSO DFT workflows.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    conv = sub.add_parser(
        "converge",
        help="Analyse ecutwfc convergence from a folder of pw.x outputs",
    )
    conv.add_argument("folder", help="Folder containing the pw.x output files")
    conv.add_argument("--pattern", default="*.out", help="Glob pattern (default: *.out)")
    conv.add_argument(
        "--threshold", type=float, default=1.0,
        help="Convergence threshold in meV/atom (default: 1.0)",
    )
    conv.add_argument("--csv", help="Write the results table to this CSV file")
    conv.add_argument("--plot", help="Save the convergence plot to this image file")

    kconv = sub.add_parser(
        "kconverge",
        help="Analyse k-point convergence from a folder of pw.x outputs (+ .in files)",
    )
    kconv.add_argument("folder", help="Folder with pw.x outputs and matching .in files")
    kconv.add_argument("--pattern", default="*.out", help="Glob pattern (default: *.out)")
    kconv.add_argument(
        "--threshold", type=float, default=1.0,
        help="Convergence threshold in meV/atom (default: 1.0)",
    )
    kconv.add_argument("--csv", help="Write the results table to this CSV file")
    kconv.add_argument("--plot", help="Save the convergence plot to this image file")
    return parser

def main(argv=None):
    args = build_parser().parse_args(argv)

    if args.command == "converge":
        try:
            points, converged = analyze_ecut(
                args.folder,
                threshold_mev=args.threshold,
                csv_path=args.csv,
                plot_path=args.plot,
                pattern=args.pattern,
            )
        except (FileNotFoundError, ValueError) as exc:
            print(f"Error: {exc}", file=sys.stderr)
            return 1

        print(f"{'file':<20}{'ecutwfc (Ry)':>14}{'E (Ry)':>18}{'dE (meV/atom)':>16}  status")
        for p in points:
            ecut = "-" if p.ecutwfc is None else f"{p.ecutwfc:g}"
            energy = "-" if p.energy_ry is None else f"{p.energy_ry:.8f}"
            delta = "-" if p.delta_mev_per_atom is None else f"{p.delta_mev_per_atom:.4f}"
            status = "ok" if p.valid else "INVALID (unfinished/unconverged)"
            print(f"{p.file:<20}{ecut:>14}{energy:>18}{delta:>16}  {status}")

        if converged is None:
            print(f"\nNot converged below {args.threshold:g} meV/atom in this range.")
        else:
            print(f"\nConverged at ecutwfc = {converged:g} Ry "
                  f"(threshold {args.threshold:g} meV/atom).")
    if args.command == "kconverge":
        try:
            points, converged = analyze_kpoints(
                args.folder,
                threshold_mev=args.threshold,
                csv_path=args.csv,
                plot_path=args.plot,
                pattern=args.pattern,
            )
        except (FileNotFoundError, ValueError) as exc:
            print(f"Error: {exc}", file=sys.stderr)
            return 1

        print(f"{'file':<20}{'k grid':>12}{'nk(IBZ)':>10}{'E (Ry)':>18}{'dE (meV/atom)':>16}  status")
        for p in points:
            ibz = "-" if p.nk_irreducible is None else str(p.nk_irreducible)
            energy = "-" if p.energy_ry is None else f"{p.energy_ry:.8f}"
            delta = "-" if p.delta_mev_per_atom is None else f"{p.delta_mev_per_atom:.4f}"
            status = "ok" if p.valid else "INVALID (unfinished/unconverged/no k grid)"
            print(f"{p.file:<20}{p.grid_label:>12}{ibz:>10}{energy:>18}{delta:>16}  {status}")

        if converged is None:
            print(f"\nNot converged below {args.threshold:g} meV/atom in this range.")
        else:
            print(f"\nConverged at k grid = {converged.grid_label} "
                  f"(threshold {args.threshold:g} meV/atom).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
