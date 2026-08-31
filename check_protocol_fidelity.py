#!/usr/bin/env python3
"""Honest protocol-fidelity gate checker (2026-06 protocol-fidelity audit).

Unlike check_simulation_health.py — which only verified that a generic UDP-over-
PointToPoint flow advanced the clock — this checker verifies that an example's
KPIs were MEASURED from a real radio stack:

  * sim_health.csv must use the honest 5-column schema
        metric,value,floor,pass,provenance
  * SINR / TBLER must carry provenance `phy-trace` (read from the SpectrumPhy),
    not a closed-form formula.
  * the radio channel must be in the packet path (channel_in_path=1).
  * transport blocks must have been decoded at the PHY (phy_rx_tb >= floor).
  * measured app throughput must clear its floor (rx_throughput_mbps).

Examples whose sim_health.csv lacks any `phy-trace` provenance row are reported
as LEGACY/COSMETIC — i.e. still on the old bolted-P2P data plane and not yet
migrated to NtnRealStackHelper. This is the CI gate for the migration.

Usage:
    python3 tools/check_protocol_fidelity.py --simTime 8 \
        --out protocol_fidelity_report.md
"""
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

# (display-name, binary, time-flag, extra-args)
# Real-stack examples register here as modules migrate (Phase 1+).
REAL_STACK_EXAMPLES = [
    ("ntn-real-stack-smoke",
     "build/contrib/ntn-traffic/examples/ns3.43-ntn-real-stack-smoke-default",
     "--simTime", "--numUes=4"),
    ("oran-ntn-real-stack-scenario",
     "build/contrib/oran-ntn/examples/ns3.43-oran-ntn-real-stack-scenario-default",
     "--duration", "--numUes=6"),
    ("ntn-cho-handover-traffic",
     "build/contrib/ntn-cho/examples/ns3.43-ntn-cho-handover-traffic-default",
     "--simSeconds", "--numUes=4"),
    ("ntn-cho-leo-basic",
     "build/contrib/ntn-cho/examples/ns3.43-ntn-cho-leo-basic-default",
     "--simTime", "--numUes=4"),
    ("ntn-cho-real-stack",
     "build/contrib/ntn-cho/examples/ns3.43-ntn-cho-real-stack-default",
     "--duration", "--numUes=4"),
    ("thz-ntn-real-stack",
     "build/contrib/thz-ntn/examples/ns3.43-thz-ntn-real-stack-default",
     "--duration", "--numUes=3"),
    ("ntn-rrc-leo-pass",
     "build/contrib/ntn-rrc/examples/ns3.43-ntn-rrc-leo-pass-default",
     "--simTime", "--numUes=4"),
    ("ntn-rrc-real-stack",
     "build/contrib/ntn-rrc/examples/ns3.43-ntn-rrc-real-stack-default",
     "--duration", "--numUes=4"),
    ("ntn-fapi-dl-data-slotloop",
     "build/contrib/ntn-fapi/examples/ns3.43-ntn-fapi-dl-data-slotloop-default",
     "--simSeconds", "--numUes=4"),
    ("ntn-fapi-real-stack",
     "build/contrib/ntn-fapi/examples/ns3.43-ntn-fapi-real-stack-default",
     "--duration", "--numUes=4"),
    ("ntn-v2x-real-stack",
     "build/contrib/ntn-v2x/examples/ns3.43-ntn-v2x-real-stack-default",
     "--duration", "--numVehicles=6"),
    ("sagin-a2g-real-stack",
     "build/contrib/ntn-sagin/examples/ns3.43-sagin-a2g-real-stack-default",
     "--duration", "--numUes=4"),
    ("ntn-sionna-cir-real-stack",
     "build/contrib/ntn-sionna/examples/ns3.43-ntn-sionna-cir-real-stack-default",
     "--duration", "--numUes=4"),
    ("ntn-slice-real-stack",
     "build/contrib/ntn-slice/examples/ns3.43-ntn-slice-real-stack-default",
     "--duration", "--numUes=9"),
    ("ntn-e2e-full-stack",
     "build/contrib/oran-ntn/examples/ns3.43-ntn-e2e-full-stack-default",
     "--duration", "--numUes=9"),
    # #197 — thz-ntn traffic examples on the channel-plugin recipe.
    ("thz-ntn-weather-traffic",
     "build/contrib/thz-ntn/examples/ns3.43-thz-ntn-weather-traffic-default",
     "--simSeconds", ""),
    ("thz-ntn-leo-ground-downlink-traffic",
     "build/contrib/thz-ntn/examples/ns3.43-thz-ntn-leo-ground-downlink-traffic-default",
     "--simSeconds", ""),
    ("thz-ntn-isl-traffic",
     "build/contrib/thz-ntn/examples/ns3.43-thz-ntn-isl-traffic-default",
     "--simSeconds", ""),
    ("thz-ntn-ris-relay-traffic",
     "build/contrib/thz-ntn/examples/ns3.43-thz-ntn-ris-relay-traffic-default",
     "--simSeconds", ""),
    ("thz-ntn-ric-controlled-traffic",
     "build/contrib/thz-ntn/examples/ns3.43-thz-ntn-ric-controlled-traffic-default",
     "--simSeconds", ""),
    ("thz-ntn-isac-coexist-traffic",
     "build/contrib/thz-ntn/examples/ns3.43-thz-ntn-isac-coexist-traffic-default",
     "--simSeconds", ""),
    ("thz-ntn-beam-tracking",
     "build/contrib/thz-ntn/examples/ns3.43-thz-ntn-beam-tracking-default",
     "--simSeconds", ""),
    # #197 — ntn-sionna traffic examples on the channel-plugin recipe.
    ("ntn-sionna-leo-downlink-traffic",
     "build/contrib/ntn-sionna/examples/ns3.43-ntn-sionna-leo-downlink-traffic-default",
     "--simSeconds", ""),
    ("ntn-sionna-rain-event-traffic",
     "build/contrib/ntn-sionna/examples/ns3.43-ntn-sionna-rain-event-traffic-default",
     "--simSeconds", ""),
    ("ntn-sionna-ris-relay-traffic",
     "build/contrib/ntn-sionna/examples/ns3.43-ntn-sionna-ris-relay-traffic-default",
     "--simSeconds", ""),
    ("ntn-sionna-mimo-traffic",
     "build/contrib/ntn-sionna/examples/ns3.43-ntn-sionna-mimo-traffic-default",
     "--simSeconds", ""),
    ("ntn-sionna-constellation-handover-traffic",
     "build/contrib/ntn-sionna/examples/ns3.43-ntn-sionna-constellation-handover-traffic-default",
     "--simSeconds", ""),
    ("ntn-sionna-composed-channel-traffic",
     "build/contrib/ntn-sionna/examples/ns3.43-ntn-sionna-composed-channel-traffic-default",
     "--simSeconds", ""),
    # #198 — slice / observability examples on the measured radio.
    ("ntn-slice-isolation-traffic",
     "build/contrib/ntn-slice/examples/ns3.43-ntn-slice-isolation-traffic-default",
     "--simSeconds", ""),
    ("ntn-three-slice-leo-geo",
     "build/contrib/ntn-slice/examples/ns3.43-ntn-three-slice-leo-geo-default",
     "--simTime", ""),
    ("ntn-observability-demo",
     "build/contrib/ntn-observability/examples/ns3.43-ntn-observability-demo-default",
     "--simTime", ""),
    # WS0 (AI-Native ORAN-NTN plan) — last SnrToPer stragglers on the real cell.
    ("oran-ntn-ric-controlled-traffic",
     "build/contrib/oran-ntn/examples/ns3.43-oran-ntn-ric-controlled-traffic-default",
     "--simSeconds", ""),
    ("ntn-observability-traffic",
     "build/contrib/ntn-observability/examples/ns3.43-ntn-observability-traffic-default",
     "--simSeconds", "--out=/tmp/ntn-obs-gate.lp"),
    # WS6 (AI-Native ORAN-NTN plan) — use-case flagships on the real cell.
    ("oran-ntn-emergency-communication",
     "build/contrib/oran-ntn/examples/ns3.43-oran-ntn-emergency-communication-default",
     "--simSeconds", "--disasterAt=3"),
    ("ntn-sagin-remote-coverage",
     "build/contrib/ntn-sagin/examples/ns3.43-ntn-sagin-remote-coverage-default",
     "--simSeconds", "--sharing=1"),
    ("ntn-v2x-edge-urllc",
     "build/contrib/ntn-v2x/examples/ns3.43-ntn-v2x-edge-urllc-default",
     "--simSeconds", "--edge=sat"),
]

# Provenance values that prove a measured (not computed) origin.
MEASURED_PROVENANCE = {"phy-trace", "packetsink", "app-trace", "flowmonitor"}


@dataclass
class Result:
    name: str
    error: str | None = None
    legacy: bool = False
    rows: dict = field(default_factory=dict)  # metric -> (value, floor, pass, provenance)
    failures: list[str] = field(default_factory=list)

    @property
    def passed(self) -> bool:
        return self.error is None and not self.legacy and not self.failures


def parse_health(path: Path) -> tuple[dict, bool]:
    """Return (rows, has_provenance). rows: metric -> (value, floor, pass, prov)."""
    rows: dict = {}
    has_prov = False
    with path.open() as f:
        reader = csv.DictReader(f)
        has_prov = "provenance" in (reader.fieldnames or [])
        for r in reader:
            prov = r.get("provenance", "")
            if prov:
                has_prov = True
            rows[r["metric"]] = (r.get("value", ""), r.get("floor", ""),
                                 r.get("pass", ""), prov)
    return rows, has_prov


def gate(rows: dict, metric: str) -> bool:
    v = rows.get(metric)
    return bool(v) and v[2] == "1"


def run_one(repo: Path, name: str, exe_rel: str, flag: str, extra: str,
            sim_time: float) -> Result:
    res = Result(name=name)
    exe = repo / exe_rel
    if not exe.exists():
        res.error = f"binary not found: {exe_rel}"
        return res
    out_dir = Path(tempfile.mkdtemp(prefix="fidelity-"))
    cmd = [str(exe), f"{flag}={sim_time}", f"--outputDir={out_dir}"]
    if extra:
        cmd += extra.split()
    try:
        subprocess.run(cmd, cwd=repo, capture_output=True, text=True,
                       timeout=max(300.0, sim_time * 30))
    except subprocess.TimeoutExpired:
        res.error = "timeout"
        shutil.rmtree(out_dir, ignore_errors=True)
        return res
    health = out_dir / "sim_health.csv"
    if not health.exists():
        res.error = "sim_health.csv missing"
        shutil.rmtree(out_dir, ignore_errors=True)
        return res
    rows, has_prov = parse_health(health)
    res.rows = rows
    shutil.rmtree(out_dir, ignore_errors=True)

    # Legacy/cosmetic detection: no phy-trace provenance anywhere.
    provs = {v[3] for v in rows.values()}
    if "phy-trace" not in provs:
        res.legacy = True
        res.failures.append("no phy-trace provenance (legacy cosmetic data plane)")
        return res

    # Honest fidelity gates.
    if not gate(rows, "channel_in_path"):
        res.failures.append("channel_in_path")
    if not gate(rows, "phy_rx_tb"):
        res.failures.append("phy_rx_tb")
    if not gate(rows, "rx_throughput_mbps"):
        res.failures.append("rx_throughput_mbps")
    # SINR/TBLER must come from a PHY trace.
    for m in ("dl_sinr_db", "dl_tbler_mean"):
        v = rows.get(m)
        if not v or v[3] != "phy-trace":
            res.failures.append(f"{m}:provenance")
    return res


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--simTime", type=float, default=8.0)
    ap.add_argument("--repo", type=Path,
                    default=Path(__file__).resolve().parent.parent)
    ap.add_argument("--out", type=Path, default=Path("protocol_fidelity_report.md"))
    ap.add_argument("--only", action="append", default=None)
    args = ap.parse_args()

    results: list[Result] = []
    for name, exe, flag, extra in REAL_STACK_EXAMPLES:
        if args.only and not any(p in name for p in args.only):
            continue
        print(f"[fidelity] running {name} ...", flush=True)
        results.append(run_one(args.repo, name, exe, flag, extra, args.simTime))

    failed = [r for r in results if not r.passed]
    with args.out.open("w") as f:
        f.write("# Protocol-Fidelity Report\n\n")
        f.write(f"- simTime: {args.simTime} s\n")
        f.write(f"- examples: {len(results)}  |  failed/legacy: {len(failed)}\n\n")
        f.write("| Example | PHY rx TB | DL SINR (dB) | TBLER | Thr (Mbps) | Provenance | Verdict |\n")
        f.write("|---|---:|---:|---:|---:|:--|:--|\n")
        for r in results:
            if r.error:
                f.write(f"| `{r.name}` | - | - | - | - | - | ERR ({r.error}) |\n")
                continue
            tb = r.rows.get("phy_rx_tb", ("-",))[0]
            sinr = r.rows.get("dl_sinr_db", ("-",))[0]
            tbler = r.rows.get("dl_tbler_mean", ("-",))[0]
            thr = r.rows.get("rx_throughput_mbps", ("-",))[0]
            prov = "phy-trace" if not r.legacy else "LEGACY"
            verdict = "REAL-STACK PASS" if r.passed else (
                "LEGACY/COSMETIC" if r.legacy else "FAIL (" + ", ".join(r.failures) + ")")
            f.write(f"| `{r.name}` | {tb} | {sinr} | {tbler} | {thr} | {prov} | {verdict} |\n")
    print(f"[fidelity] report -> {args.out}")
    if failed:
        for r in failed:
            print(f"[fidelity] FAIL {r.name}: {r.error or ', '.join(r.failures)}")
        return 1
    print("[fidelity] PASS — all checked examples are REAL-STACK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
