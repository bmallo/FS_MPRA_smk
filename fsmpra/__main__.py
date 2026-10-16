"""Run with python -m fsmpra from the repository root."""

import argparse
import json
import logging


def main():
    parser = argparse.ArgumentParser(description="LDLR Fiber-seq redesign pilot (descriptive measurements)")
    commands = parser.add_subparsers(dest="command", required=True)
    catalog = commands.add_parser("catalog", help="Validate design, coordinates, inputs, and published MPRA joins")
    catalog.add_argument("--manifest", required=True)
    catalog.add_argument("--output", required=True)
    molecules = commands.add_parser("molecules", help="Sample the complete ledger and join source metrics")
    molecules.add_argument("--catalog", required=True)
    molecules.add_argument("--sample", required=True)
    molecules.add_argument("--output", required=True)
    molecules.add_argument("--fraction", type=float, default=0.01)
    molecules.add_argument("--seed", type=int, default=42)
    molecules.add_argument("--max-molecules", type=int, default=200000)
    features = commands.add_parser("measure", help="Extract small-footprint and MSP accessibility measurements")
    features.add_argument("--catalog", required=True)
    features.add_argument("--molecules", required=True)
    features.add_argument("--output", required=True)
    features.add_argument("--min-passes", type=int, default=10)
    features.add_argument("--tf-min", type=int, default=10)
    features.add_argument("--tf-max", type=int, default=60)
    features.add_argument("--nuc-min", type=int, default=85)
    features.add_argument("--mask-radius", type=int, default=2)
    features.add_argument("--examples-per-variant", type=int, default=4)
    plot = commands.add_parser("plot", help="Export a descriptive scientific overview (requires matplotlib and numpy)")
    plot.add_argument("--pilot", required=True)
    plot.add_argument("--output", required=True)
    plot.add_argument("--min-reads", type=int, default=10)
    dashboard = commands.add_parser("dashboard", help="Export a local self-contained interactive dashboard")
    dashboard.add_argument("--pilot", required=True)
    dashboard.add_argument("--output", required=True)
    comparison = commands.add_parser("compare", help="Compare runs with explicit reaction conditions; no pooling")
    comparison.add_argument("--pilots", nargs="+", required=True)
    comparison.add_argument("--manifest", required=True)
    comparison.add_argument("--output", required=True)
    comparison.add_argument("--min-reads", type=int, default=10)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    if args.command == "catalog":
        from .catalog import build_catalog
        result = build_catalog(args.manifest, args.output)
        result = {"libraries": [{"id": l["id"], "class_counts": l["class_counts"]} for l in result["libraries"]],
                  "mpra_audit": result["mpra_audit"], "output": args.output}
    elif args.command == "molecules":
        from .molecules import import_molecules
        result = import_molecules(args.catalog, args.sample, args.output, args.fraction, args.seed, args.max_molecules)
    elif args.command == "measure":
        from .features import measure
        result = measure(args.catalog, args.molecules, args.output, args.min_passes, args.tf_min,
                         args.tf_max, args.nuc_min, args.mask_radius, args.examples_per_variant)
    elif args.command == "plot":
        from .plot import plot_pilot
        result = plot_pilot(args.pilot, args.output, args.min_reads)
    elif args.command == "compare":
        from .compare import compare_pilots
        result = compare_pilots(args.pilots, args.manifest, args.output, args.min_reads)
    else:
        from .dashboard import export_dashboard
        result = export_dashboard(args.pilot, args.output)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
