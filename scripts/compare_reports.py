#!/usr/bin/env python3
"""
Compare two pipeline runs and generate a comparative report (HTML + plots).
"""

from __future__ import annotations

import argparse
from pathlib import Path

from orchestrator.visualization.reports import generate_comparative_report


def main():
    p = argparse.ArgumentParser(description="Generate comparative report for two runs")
    p.add_argument("--run1-results", required=True, help="Path to results.csv for run A")
    p.add_argument("--run1-manifest", required=False, help="Path to manifest.json for run A")
    p.add_argument("--run2-results", required=True, help="Path to results.csv for run B")
    p.add_argument("--run2-manifest", required=False, help="Path to manifest.json for run B")
    p.add_argument("--out", required=False, help="Output directory", default="data/outputs/compare")
    p.add_argument("--title", required=False, help="Report title", default="Comparative Report")
    args = p.parse_args()

    out_dir = Path(args.out)
    files = generate_comparative_report(
        Path(args.run1_results),
        Path(args.run1_manifest) if args.run1_manifest else None,
        Path(args.run2_results),
        Path(args.run2_manifest) if args.run2_manifest else None,
        out_dir,
        title=args.title,
    )
    print(f"Comparative report: {files.get('html')}")


if __name__ == "__main__":
    main()

