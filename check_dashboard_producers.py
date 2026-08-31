#!/usr/bin/env python3
"""Every measurement a shipped dashboard queries must have a real producer.

OBS-03: two of the four Grafana dashboards queried measurements that nothing in
the toolkit ever wrote. `ntn_isl` and `ntn_handover` were declared in
`ntn-metric-schema.h` and referenced nowhere else except a unit test that
asserted the string constant equalled its own literal, which is a test that
cannot fail. Both dashboards therefore rendered empty forever, and nothing in
the build said so.

This checks the chain end to end:

  dashboard JSON  ->  declared in ntn-metric-schema.h  ->  pushed by real code

A measurement referenced only from `test/` counts as UNPRODUCED, because that is
exactly the state the defect was in.

Exit status is 0 when every queried measurement has a producer, 1 otherwise.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import re
import sys

# Measurement names look like ntn_<lowercase_with_underscores>.
MEASUREMENT_RE = re.compile(r"\bntn_[a-z][a-z0-9_]*\b")


def dashboard_measurements(dash_dir: pathlib.Path) -> dict[str, set[str]]:
    """Measurement names referenced by each dashboard file."""
    out: dict[str, set[str]] = {}
    for path in sorted(dash_dir.rglob("*.json")):
        try:
            text = path.read_text(encoding="utf-8", errors="replace")
            json.loads(text)  # reject anything that is not really a dashboard
        except (OSError, json.JSONDecodeError):
            continue
        names = set(MEASUREMENT_RE.findall(text))
        if names:
            out[path.name] = names
    return out


def schema_names(schema: pathlib.Path) -> dict[str, str]:
    """Map measurement string -> C++ constant name, from the schema header."""
    text = schema.read_text(encoding="utf-8", errors="replace")
    pairs = re.findall(
        r'inline\s+constexpr\s+const\s+char\*\s+(\w+)\s*=\s*"(ntn_[a-z0-9_]+)"', text
    )
    return {value: const for const, value in pairs}


def producer_sites(contrib: pathlib.Path, const: str) -> list[str]:
    """Non-test source lines that push this measurement."""
    hits: list[str] = []
    needle = f"measurement::{const}"
    for path in contrib.rglob("*.cc"):
        parts = path.parts
        if "test" in parts:
            continue  # a test asserting the constant is not a producer
        if path.name.endswith("-test-suite.cc"):
            continue
        try:
            text = path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        if needle in text:
            hits.append(str(path))
    return hits


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", default=".", help="ns-3-dev root")
    ap.add_argument(
        "--dashboards",
        default="../ntn-control-center/docker/grafana-provisioning/dashboards",
        help="directory of shipped Grafana dashboards",
    )
    args = ap.parse_args()

    repo = pathlib.Path(args.repo).resolve()
    dash_dir = (repo / args.dashboards).resolve()
    schema = repo / "contrib/ntn-observability/model/ntn-metric-schema.h"
    contrib = repo / "contrib"

    if not schema.is_file():
        print(f"[dash] schema not found: {schema}", file=sys.stderr)
        return 1
    if not dash_dir.is_dir():
        print(f"[dash] no dashboards at {dash_dir}; nothing to check")
        return 0

    known = schema_names(schema)
    dashboards = dashboard_measurements(dash_dir)
    if not dashboards:
        print(f"[dash] no dashboard JSON under {dash_dir}; nothing to check")
        return 0

    undeclared: list[tuple[str, str]] = []
    unproduced: list[tuple[str, str]] = []
    ok = 0

    for dash, names in dashboards.items():
        for name in sorted(names):
            const = known.get(name)
            if const is None:
                undeclared.append((dash, name))
                continue
            if not producer_sites(contrib, const):
                unproduced.append((dash, name))
            else:
                ok += 1

    print(f"[dash] dashboards: {len(dashboards)}  measurements OK: {ok}")
    for dash, name in undeclared:
        print(f"[dash] ERROR {dash}: '{name}' is not declared in ntn-metric-schema.h")
    for dash, name in unproduced:
        print(
            f"[dash] ERROR {dash}: '{name}' is declared but NOTHING pushes it "
            f"(outside test/), so this panel renders empty"
        )

    if undeclared or unproduced:
        print("[dash] FAIL")
        return 1
    print("[dash] PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
