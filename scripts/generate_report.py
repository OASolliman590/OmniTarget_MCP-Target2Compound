#!/usr/bin/env python3
"""
Generate a static HTML report from a pipeline run directory or explicit files.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from orchestrator.visualization.reports import generate_report


def main():
    p = argparse.ArgumentParser(description="Generate pipeline report")
    p.add_argument("--results", required=True, help="Path to results.csv")
    p.add_argument("--manifest", required=False, help="Path to manifest.json")
    p.add_argument("--out", required=False, help="Output report directory", default=None)
    p.add_argument("--title", required=False, help="Report title", default="Pipeline Report")
    args = p.parse_args()

    results_csv = Path(args.results)
    manifest = Path(args.manifest) if args.manifest else None
    out_dir = Path(args.out) if args.out else results_csv.parent / "report"
    files = generate_report(results_csv, manifest, out_dir, title=args.title)
    print(f"Report generated: {files.get('html')}\n")


if __name__ == "__main__":
    main()

