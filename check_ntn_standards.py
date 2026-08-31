#!/usr/bin/env python3
# Copyright (c) 2026 Muhammad Uzair
# SPDX-License-Identifier: GPL-2.0-only
"""WS5 standards-validation gate (AI-Native ORAN-NTN adoption).

Aggregates every standards check into one PASS/FAIL:
  1. unit/e2e test suites (orbital theory, ORAN apps, KPM monitor, multi-tier
     RIC, WS4 payload/split/role-switch) via ./test.py;
  2. TR 38.821 Set-1 LEO-600 CNR calibration on the measured radio
     (offset constant across the pass + elevation slope matches FSPL);
  3. paper-Table-3 per-platform latency bands with real packets;
  4. the five 3GPP NTN handover trigger classes all execute a handover on
     the measured radio (a3 / d1 / t1 / elevation / ta).

Run from the ns-3 root:  python3 tools/check_ntn_standards.py
"""

import re
import subprocess
import sys

NS3_ROOT = "."

TEST_SUITES = [
    "ntn-standards-validation",
    "ntn-oran-application",
    "ntn-oran-ai-flow-monitor",
    "oran-ntn-multi-tier-ric",
    "oran-ntn-ws4",
]

failures = []


def run(cmd, timeout=900):
    return subprocess.run(cmd, shell=True, cwd=NS3_ROOT, capture_output=True,
                          text=True, timeout=timeout)


def check(name, ok, detail=""):
    print(f"[standards] {'PASS' if ok else 'FAIL'}  {name}  {detail}")
    if not ok:
        failures.append(name)


def main():
    # 1. test suites
    for suite in TEST_SUITES:
        r = run(f"./test.py -s {suite}")
        ok = "FAIL" not in r.stdout and "CRASH" not in r.stdout and r.returncode == 0
        check(f"test-suite {suite}", ok)

    # 2. TR 38.821 calibration (exit code carries the verdict)
    r = run('./ns3 run "ntn-tr38821-calibration --simSeconds=120"')
    # WF-08 gate 5. The gate is stated as "calibration offset ~ array gain,
    # offset_std < 1.5 dB". This used to parse offset_std and then assert only
    # that the word CALIBRATED appeared in the output - so the numeric bound the
    # gate is named for was never checked, and a run with a 5 dB spread would
    # have passed as long as the example still printed its verdict string.
    m = re.search(r"offset_std=([\d.]+) dB", r.stdout)
    m_mean = re.search(r"offset_mean=([-\d.]+) dB", r.stdout)
    m_delta = re.search(r"meas_delta=([-\d.]+)\s+tr_delta=([-\d.]+)\s+\(tol=([\d.]+)\)",
                        r.stdout)
    std = float(m.group(1)) if m else float("inf")
    mean = float(m_mean.group(1)) if m_mean else float("nan")
    reasons = []
    if r.returncode != 0:
        reasons.append("nonzero exit")
    if "CALIBRATED" not in r.stdout:
        reasons.append("no CALIBRATED verdict")
    if not m:
        reasons.append("offset_std not reported")
    elif std >= 1.5:
        reasons.append(f"offset_std {std:.2f} dB >= 1.5 dB bound")
    if m_delta:
        # The elevation SLOPE must track TR 38.821 within the example's own
        # stated tolerance: a constant offset can be calibrated away, a wrong
        # slope cannot.
        meas_d, tr_d, tol = (float(m_delta.group(1)), float(m_delta.group(2)),
                             float(m_delta.group(3)))
        if abs(meas_d - tr_d) > tol:
            reasons.append(
                f"elevation slope deviates {abs(meas_d - tr_d):.2f} dB > tol {tol:.2f} dB")
    else:
        reasons.append("elevation slope not reported")
    check("tr38821-calibration", not reasons,
          f"offset_mean={mean:.2f} dB offset_std={std:.2f} dB"
          + (" | " + "; ".join(reasons) if reasons else ""))

    # CVC-13: user-facing capability claims must not contradict the code.
    r = run("python3 tools/check_doc_claims.py")
    check("doc-claims", r.returncode == 0,
          (r.stdout.strip().splitlines() or ["no output"])[0])

    # The audit register's status tags must agree with the implementation log
    # and with the fix markers in the tree. Stale tags are not cosmetic: they
    # sent more than one session off re-deriving closures that were already
    # recorded, so the reconciliation runs as a gate rather than by hand.
    r = run("python3 tools/check_gap_register.py --quiet")
    note = (r.stdout.strip().splitlines() or ["no output"])[-1].replace("[gap-register] ", "")
    if note.startswith("SKIP"):
        # The engineering records this gate reconciles are kept out of the public
        # tree, so a clone will not have them. Say skipped rather than passed: a
        # gate that reports PASS when it never ran is the exact habit the gate
        # exists to catch.
        print(f"[standards] SKIP  gap-register  {note}")
    else:
        check("gap-register", r.returncode == 0, note)

    # 3. Table-3 latency bands
    r = run("./ns3 run ntn-platform-latency-validation")
    check("table3-latency-bands", r.returncode == 0 and "ALL CLASSES PASS" in r.stdout)

    # 4. NTN trigger classes: Rel-17 CondEvents A4-combined a3/d1/t1,
    #    Rel-18 CondEventD2 (moving references), plus the TR 38.821-studied
    #    elevation and timing-advance mechanisms.
    for trig in ["a3", "d1", "t1", "d2", "elevation", "ta"]:
        r = run(f'./ns3 run "ntn-cho-handover-traffic --simSeconds=60 --numUes=2 '
                f'--trigger={trig}"')
        m = re.search(r"handovers=(\d+)", r.stdout)
        handovers = int(m.group(1)) if m else 0
        check(f"ho-trigger-class {trig}", r.returncode == 0 and handovers >= 1,
              f"handovers={handovers}")

    print()
    if failures:
        print(f"[standards] FAIL — {len(failures)} check(s): {', '.join(failures)}")
        return 1
    print("[standards] PASS — toolkit validated against 3GPP TR 38.821 link "
          "budget, orbital theory, paper Table 3, and all five NTN HO classes")
    return 0


if __name__ == "__main__":
    sys.exit(main())
