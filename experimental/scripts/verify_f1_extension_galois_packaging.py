#!/usr/bin/env python3
"""Galois-package the F1 extension full-orbit scan signal.

Consumes the responsible F_7^4/F_7 pole-conjugate-pair row from
``verify_f1_extension_full_orbit_scan.py`` and packages the Frobenius orbit of
the fixed received line.  The question is whether the raw K=F signal collapses
to base/tower-paid strata or becomes an orbit-closed K=F inventory candidate.
"""

from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path
from typing import Any

import verify_f1_extension_full_orbit_scan as raw


REPO = Path(__file__).resolve().parents[2]
ARTIFACT = (
    REPO
    / "experimental"
    / "data"
    / "certificates"
    / "f1-extension-galois-packaging"
    / "f1_extension_galois_packaging.json"
)


def conjugate_values(
    field: raw.Field,
    values: dict[tuple[int, ...], tuple[int, ...]],
    steps: int,
) -> dict[tuple[int, ...], tuple[int, ...]]:
    return {x: field.frobenius(value, steps) for x, value in values.items()}


def degree_buckets(field: raw.Field, slopes: set[tuple[int, ...]]) -> dict[str, int]:
    return {
        str(d): sum(1 for slope in slopes if raw.minimal_field_degree(field, slope) == d)
        for d in raw.divisors(field.degree)
    }


def bad_slope_set(
    field: raw.Field,
    domain: list[tuple[int, ...]],
    k: int,
    agreement: int,
    f_values: dict[tuple[int, ...], tuple[int, ...]],
    g_values: dict[tuple[int, ...], tuple[int, ...]],
) -> set[tuple[int, ...]]:
    supports = [
        tuple(support)
        for size in range(agreement, len(domain) + 1)
        for support in combinations(domain, size)
    ]
    f_good = {
        support: raw.degree_lt_k_on_support(field, domain, f_values, support, k)
        for support in supports
    }
    g_good = {
        support: raw.degree_lt_k_on_support(field, domain, g_values, support, k)
        for support in supports
    }
    bad: set[tuple[int, ...]] = set()
    for slope in field.elements():
        h_values = {
            x: field.add(f_values[x], field.mul(slope, g_values[x]))
            for x in domain
        }
        for support in supports:
            if f_good[support] and g_good[support]:
                continue
            if raw.degree_lt_k_on_support(field, domain, h_values, support, k):
                bad.add(slope)
                break
    return bad


def slope_payload(field: raw.Field, slopes: set[tuple[int, ...]]) -> dict[str, Any]:
    full = sorted(
        slope for slope in slopes if raw.minimal_field_degree(field, slope) == field.degree
    )
    return {
        "total_bad_slopes": len(slopes),
        "minimal_field_degree_counts": degree_buckets(field, slopes),
        "full_degree_bad_slopes": len(full),
        "frobenius_orbit_diagnostic": raw.full_orbit_diagnostic(field, full),
    }


def scan_alpha(alpha_index: int) -> dict[str, Any]:
    p, extension_degree = 7, 4
    field = raw.Field(p, extension_degree, raw.find_irreducible(p, extension_degree))
    domain = [field.embed(x) for x in range(1, p)]
    n, k, agreement = p - 1, (p - 1) // 2, 4
    alphas = [
        x
        for x in field.elements()
        if x != field.zero and raw.minimal_field_degree(field, x) == extension_degree
    ]
    alpha = alphas[alpha_index]
    f0, g0 = raw.candidate_values(field, "pole_conjugate_pair", alpha, domain)

    component_sets = []
    component_reports = []
    for i in range(extension_degree):
        f_i = conjugate_values(field, f0, i)
        g_i = conjugate_values(field, g0, i)
        slopes_i = bad_slope_set(field, domain, k, agreement, f_i, g_i)
        component_sets.append(slopes_i)
        expected = {field.frobenius(slope, i) for slope in component_sets[0]}
        report = {
            "component": i,
            "equals_frobenius_conjugate_of_component_0": slopes_i == expected,
            **slope_payload(field, slopes_i),
        }
        component_reports.append(report)

    union_slopes = set().union(*component_sets)
    intersection_slopes = set.intersection(*component_sets)
    union_payload = slope_payload(field, union_slopes)
    full_orbit_closed = union_payload["frobenius_orbit_diagnostic"][
        "full_orbit_closed_for_fixed_line"
    ]
    full_orbits = union_payload["frobenius_orbit_diagnostic"]["complete_full_orbits"]
    return {
        "case": f"F7_e4_alpha{alpha_index}_pole_conjugate_pair_galois_package",
        "p": p,
        "extension_degree": extension_degree,
        "n": n,
        "k": k,
        "agreement": agreement,
        "alpha": field.encode(alpha),
        "component_reports": component_reports,
        "union": union_payload,
        "intersection": slope_payload(field, intersection_slopes),
        "classification": {
            "collapses_to_base_or_tower_paid_strata": union_payload[
                "full_degree_bad_slopes"
            ]
            == 0,
            "orbit_closed_K_eq_F_inventory_candidate": (
                union_payload["full_degree_bad_slopes"] > 0 and full_orbit_closed
            ),
            "complete_K_eq_F_orbits": full_orbits,
            "raw_fixed_line_signal_was_not_orbit_closed": not component_reports[0][
                "frobenius_orbit_diagnostic"
            ]["full_orbit_closed_for_fixed_line"],
            "inventory_reading": (
                "The Galois package does not collapse to paid base/tower strata: "
                "the union carries 14 complete full-degree Frobenius orbits "
                "(56 K=F slopes). This is an orbit-closed K=F inventory "
                "candidate for the extension-cell branch, still at toy scale."
            ),
        },
    }


def build_certificate() -> dict[str, Any]:
    cases = [scan_alpha(0), scan_alpha(1)]
    cert = {
        "schema": "f1-extension-galois-packaging-v1",
        "status": "EXPERIMENTAL_AUDIT",
        "source_packet": "f1-extension-full-orbit-scan-v1",
        "target": (
            "Galois-package the F_7^4/F_7 pole-conjugate-pair signal and "
            "decide whether it collapses to base/tower-paid strata."
        ),
        "cases": cases,
        "summary": {
            "cases_checked": len(cases),
            "all_components_are_frobenius_conjugates": all(
                report["equals_frobenius_conjugate_of_component_0"]
                for case in cases
                for report in case["component_reports"]
            ),
            "max_union_full_degree_bad_slopes": max(
                case["union"]["full_degree_bad_slopes"] for case in cases
            ),
            "max_complete_K_eq_F_orbits": max(
                case["classification"]["complete_K_eq_F_orbits"] for case in cases
            ),
            "any_case_collapses_to_base_or_tower": any(
                case["classification"]["collapses_to_base_or_tower_paid_strata"]
                for case in cases
            ),
            "all_cases_are_orbit_closed_K_eq_F_inventory_candidates": all(
                case["classification"]["orbit_closed_K_eq_F_inventory_candidate"]
                for case in cases
            ),
            "reading": (
                "The responsible raw fixed-line signal becomes an orbit-closed "
                "K=F inventory candidate after Galois packaging: 56 full-degree "
                "slopes, arranged as 14 complete Frobenius orbits, for both "
                "checked alpha choices. This is toy-scale evidence for the "
                "extension-cell blocker, not a deployed counterexample."
            ),
        },
        "non_claims": [
            "does not prove asymptotic growth of the K=F inventory",
            "does not classify all F-valued pairs",
            "does not prove the descended-cycle inventory theorem",
            "does not pay or refute the deployed paid_extension cell",
        ],
    }
    validate(cert)
    return cert


def validate(cert: dict[str, Any]) -> None:
    if cert.get("schema") != "f1-extension-galois-packaging-v1":
        raise AssertionError("unexpected schema")
    if cert.get("status") != "EXPERIMENTAL_AUDIT":
        raise AssertionError("unexpected status")
    if len(cert.get("cases", [])) != 2:
        raise AssertionError("expected two alpha cases")
    for case in cert["cases"]:
        if case["union"]["total_bad_slopes"] != 58:
            raise AssertionError(f"unexpected union size in {case['case']}")
        if case["union"]["full_degree_bad_slopes"] != 56:
            raise AssertionError(f"unexpected K=F union size in {case['case']}")
        if case["union"]["minimal_field_degree_counts"] != {"1": 0, "2": 2, "4": 56}:
            raise AssertionError(f"unexpected union buckets in {case['case']}")
        diag = case["union"]["frobenius_orbit_diagnostic"]
        if not diag["full_orbit_closed_for_fixed_line"]:
            raise AssertionError(f"union not orbit-closed in {case['case']}")
        if diag["complete_full_orbits"] != 14 or diag["touched_full_orbits"] != 14:
            raise AssertionError(f"unexpected orbit count in {case['case']}")
        if case["intersection"]["total_bad_slopes"] != 0:
            raise AssertionError(f"unexpected common-slope intersection in {case['case']}")
        if not all(
            report["equals_frobenius_conjugate_of_component_0"]
            for report in case["component_reports"]
        ):
            raise AssertionError(f"component conjugacy failed in {case['case']}")
    if not cert["summary"]["all_cases_are_orbit_closed_K_eq_F_inventory_candidates"]:
        raise AssertionError("summary lost inventory-candidate status")


def assert_same(expected: dict[str, Any], actual: dict[str, Any]) -> None:
    if expected != actual:
        raise AssertionError(
            "certificate mismatch\nexpected:\n"
            + json.dumps(expected, indent=2, sort_keys=True)
            + "\nactual:\n"
            + json.dumps(actual, indent=2, sort_keys=True)
        )


def tamper_selftest(cert: dict[str, Any]) -> None:
    bad = json.loads(json.dumps(cert))
    bad["cases"][0]["union"]["full_degree_bad_slopes"] -= 1
    try:
        validate(bad)
    except AssertionError:
        return
    raise AssertionError("tamper selftest failed")


def print_summary(cert: dict[str, Any]) -> None:
    print("f1_extension_galois_packaging")
    print(f"  status: {cert['status']}")
    print(f"  reading: {cert['summary']['reading']}")
    for case in cert["cases"]:
        print(
            f"  {case['case']}: union={case['union']['total_bad_slopes']} "
            f"K=F={case['union']['full_degree_bad_slopes']} "
            f"orbits={case['classification']['complete_K_eq_F_orbits']} "
            f"intersection={case['intersection']['total_bad_slopes']}"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--emit-defaults", action="store_true", help="write the frozen certificate")
    parser.add_argument("--check", type=Path, help="check an existing certificate")
    parser.add_argument("--json", action="store_true", help="print JSON to stdout")
    parser.add_argument("--tamper-selftest", action="store_true", help="verify validator rejects a mutation")
    args = parser.parse_args()

    cert = build_certificate()
    if args.tamper_selftest:
        tamper_selftest(cert)
    if args.emit_defaults:
        ARTIFACT.parent.mkdir(parents=True, exist_ok=True)
        ARTIFACT.write_text(json.dumps(cert, indent=2, sort_keys=True) + "\n")
        print(f"wrote {ARTIFACT.relative_to(REPO)}")
    if args.check:
        actual = json.loads(args.check.read_text())
        validate(actual)
        assert_same(cert, actual)
        print(f"checked {args.check}")
    if args.json:
        print(json.dumps(cert, indent=2, sort_keys=True))
    if not args.emit_defaults and not args.check and not args.json:
        print_summary(cert)


if __name__ == "__main__":
    main()
