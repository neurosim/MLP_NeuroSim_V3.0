#!/usr/bin/env python3
"""Walks every example known to the realistic-traffic helper, runs it for
a short simulation time, and aggregates sim_health.csv into a Markdown
report.  Fails (exit 1) if any example misses its realism gates.

Usage:
    python3 tools/check_simulation_health.py --simTime 30 \
        --build-dir build --out report.md
"""
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

# (cli-name, example-binary, time-flag, default-args)
EXAMPLES = [
    ("ntn-cho-leo-basic",         "build/contrib/ntn-cho/examples/ns3.43-ntn-cho-leo-basic-default",         "--simTime",     "--numUes=8"),
    ("ntn-cho-full-constellation","build/contrib/ntn-cho/examples/ns3.43-ntn-cho-full-constellation-default","--simTime",     "--numUes=15 --algorithm=tte-aware"),
    ("oran-ntn-full-scenario",    "build/contrib/oran-ntn/examples/ns3.43-oran-ntn-full-scenario-default",    "--duration",    "--numUes=10"),
    ("thz-ntn-leo-ground",        "build/contrib/thz-ntn/examples/ns3.43-thz-ntn-leo-ground-default",         "--simTime",     ""),
    ("thz-ntn-isl",               "build/contrib/thz-ntn/examples/ns3.43-thz-ntn-isl-default",                "--simTime",     ""),
    ("thz-ntn-beam-tracking",     "build/contrib/thz-ntn/examples/ns3.43-thz-ntn-beam-tracking-default",      "--passDuration", ""),
    ("ntn-rrc-leo-pass",          "build/contrib/ntn-rrc/examples/ns3.43-ntn-rrc-leo-pass-default",           "--simTime",     ""),
    ("ntn-rrc-full-stack",        "build/contrib/ntn-rrc/examples/ns3.43-ntn-rrc-full-stack-default",         "--simTime",     ""),
    ("sagin-haps-leo-relay",      "build/contrib/ntn-sagin/examples/ns3.43-sagin-haps-leo-relay-default",     "--simTime",     ""),
    ("sagin-uav-swarm",           "build/contrib/ntn-sagin/examples/ns3.43-sagin-uav-swarm-default",          "--simTime",     ""),
    ("sagin-aeronautical",        "build/contrib/ntn-sagin/examples/ns3.43-sagin-aeronautical-default",       "--simTime",     ""),
    ("ntn-three-slice-leo-geo",   "build/contrib/ntn-slice/examples/ns3.43-ntn-three-slice-leo-geo-default",  "--simTime",     ""),
    ("ntn-v2x-rural-highway",     "build/contrib/ntn-v2x/examples/ns3.43-ntn-v2x-rural-highway-default",       "--simTime",     ""),
    ("ntn-observability-demo",    "build/contrib/ntn-observability/examples/ns3.43-ntn-observability-demo-default", "--simTime", ""),
]


@dataclass
class Health:
    name: str
    sim_sec: float = 0.0
    wall_sec: float = 0.0
    tx: int = 0
    rx: int = 0
    tx_per_sim_sec: float = 0.0
    rx_over_tx: float = 0.0
    passed: bool = False
    failures: list[str] = None
    error: str | None = None

    def __post_init__(self):
        if self.failures is None:
            self.failures = []


def parse_health(csv_path: Path) -> dict[str, tuple[float, float, int]]:
    """Returns metric -> (value, floor, pass)."""
    out: dict[str, tuple[float, float, int]] = {}
    with csv_path.open() as f:
        for row in csv.DictReader(f):
            try:
                val = float(row["value"]) if row["value"].replace(".", "", 1).replace("-", "", 1).isdigit() else row["value"]
            except Exception:
                val = row["value"]
            try:
                floor = float(row["floor"]) if row["floor"] not in ("-", "") else None
            except Exception:
                floor = None
            try:
                ok = int(row["pass"])
            except Exception:
                ok = 1
            out[row["metric"]] = (val, floor, ok)
    return out


def run_one(repo: Path, exe_rel: str, time_flag: str, extra: str,
            sim_time: float) -> Health:
    name = Path(exe_rel).name.replace("ns3.43-", "").replace("-default", "")
    h = Health(name=name)
    exe = repo / exe_rel
    if not exe.exists():
        h.error = f"binary not found: {exe}"
        return h
    out_dir = Path(tempfile.mkdtemp(prefix="simhealth-"))
    cmd = [str(exe), f"{time_flag}={sim_time}", f"--outputDir={out_dir}"]
    if extra:
        cmd += extra.split()
    proc = subprocess.run(cmd, cwd=repo, capture_output=True, text=True,
                          timeout=max(120.0, sim_time * 4))
    health_csv = out_dir / "sim_health.csv"
    if not health_csv.exists():
        # Some examples don't take outputDir for the run; try cwd as default
        for guess in [
            Path("./ntn-rrc-leo-pass-out"),
            Path(f"./{name}-out"),
        ]:
            cand = repo / guess / "sim_health.csv"
            if cand.exists():
                health_csv = cand
                break
    if not health_csv.exists():
        h.error = "sim_health.csv missing"
        h.failures.append("sim_health.csv missing")
        return h
    metrics = parse_health(health_csv)
    h.sim_sec = float(metrics["sim_time_s"][0])
    h.wall_sec = float(metrics["wall_clock_s"][0])
    h.tx = int(metrics["packets_tx"][0])
    h.rx = int(metrics["packets_rx"][0])
    h.tx_per_sim_sec = float(metrics["tx_per_sim_sec"][0])
    h.rx_over_tx = float(metrics["rx_over_tx"][0])
    for metric in ("wall_clock_s", "packets_tx", "tx_per_sim_sec", "rx_over_tx"):
        if metric in metrics and metrics[metric][2] != 1:
            h.failures.append(metric)
    h.passed = not h.failures
    shutil.rmtree(out_dir, ignore_errors=True)
    return h


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--simTime", type=float, default=30.0,
                    help="Simulation duration to pass to each example")
    ap.add_argument("--repo", type=Path, default=Path(__file__).resolve().parent.parent,
                    help="ns-3-dev root")
    ap.add_argument("--out", type=Path, default=Path("sim_health_report.md"))
    ap.add_argument("--only", action="append", default=None,
                    help="Run only examples whose name contains this substring")
    args = ap.parse_args()

    rows: list[Health] = []
    for name, exe, flag, extra in EXAMPLES:
        if args.only and not any(p in name for p in args.only):
            continue
        print(f"[check] running {name} ...", flush=True)
        rows.append(run_one(args.repo, exe, flag, extra, args.simTime))

    failed = [r for r in rows if not r.passed and r.error is None]
    skipped = [r for r in rows if r.error]
    summary_path = args.out
    with summary_path.open("w") as f:
        f.write("# Simulation Health Report\n\n")
        f.write(f"- simTime: {args.simTime} s\n")
        f.write(f"- examples run: {len(rows)}\n")
        f.write(f"- failed: {len(failed)}\n")
        f.write(f"- skipped: {len(skipped)}\n\n")
        f.write("| Example | Sim s | Wall s | Tx | Rx/Tx | Tx/SimSec | Gates |\n")
        f.write("|---|---:|---:|---:|---:|---:|:---:|\n")
        for r in rows:
            if r.error:
                f.write(f"| `{r.name}` | - | - | - | - | - | ERR ({r.error}) |\n")
                continue
            status = "PASS" if r.passed else "FAIL (" + ", ".join(r.failures) + ")"
            f.write(f"| `{r.name}` | {r.sim_sec:.0f} | {r.wall_sec:.2f} | "
                    f"{r.tx} | {r.rx_over_tx:.3f} | {r.tx_per_sim_sec:.1f} | {status} |\n")
    print(f"[check] report written to {summary_path}")
    if failed or skipped:
        print(f"[check] FAIL: {len(failed)} failed, {len(skipped)} skipped")
        return 1
    print("[check] PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
