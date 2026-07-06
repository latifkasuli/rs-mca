#!/usr/bin/env python3
"""Quartic growth scan for the F1 extension Galois-packaged signal.

This specializes the full-orbit/Galois-packaging scan to the quartic
``pole_conjugate_pair`` family and uses support-linear constraints to avoid a
full slope-by-support sweep.  It tests whether the orbit-closed K=F inventory
candidate found at F_7^4/F_7 stays bounded or grows with the base prime.
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
    / "f1-extension-quartic-growth"
    / "f1_extension_quartic_growth.json"
)
FRONTIER_NUMERATOR = 1_116_048
FRONTIER_DENOMINATOR = 2_097_152


def frontier_agreement(n: int, k: int) -> int:
    return max(
        k + 1,
        (FRONTIER_NUMERATOR * n + FRONTIER_DENOMINATOR - 1)
        // FRONTIER_DENOMINATOR,
    )


def lagrange_coefficients(
    field: raw.Field,
    base: tuple[tuple[int, ...], ...],
    x: tuple[int, ...],
) -> list[tuple[int, ...]]:
    coeffs = []
    for i, xi in enumerate(base):
        num = field.one
        den = field.one
        for j, xj in enumerate(base):
            if i == j:
                continue
            num = field.mul(num, field.sub(x, xj))
            den = field.mul(den, field.sub(xi, xj))
        coeffs.append(field.div(num, den))
    return coeffs


def linear_predict(
    field: raw.Field,
    values: dict[tuple[int, ...], tuple[int, ...]],
    base: tuple[tuple[int, ...], ...],
    coeffs: list[tuple[int, ...]],
) -> tuple[int, ...]:
    out = field.zero
    for coeff, point in zip(coeffs, base):
        out = field.add(out, field.mul(coeff, values[point]))
    return out


def support_good_and_slope(
    field: raw.Field,
    f_values: dict[tuple[int, ...], tuple[int, ...]],
    g_values: dict[tuple[int, ...], tuple[int, ...]],
    support: tuple[tuple[int, ...], ...],
    k: int,
) -> tuple[bool, bool, str, tuple[int, ...] | None]:
    base = support[:k]
    f_good = True
    g_good = True
    candidate: tuple[int, ...] | None = None
    all_slope = True
    for x in support[k:]:
        coeffs = lagrange_coefficients(field, base, x)
        f_err = field.sub(linear_predict(field, f_values, base, coeffs), f_values[x])
        g_err = field.sub(linear_predict(field, g_values, base, coeffs), g_values[x])
        if f_err != field.zero:
            f_good = False
        if g_err != field.zero:
            g_good = False
        if g_err == field.zero:
            if f_err != field.zero:
                return f_good, g_good, "none", None
            continue
        all_slope = False
        z = field.div(field.neg(f_err), g_err)
        if candidate is None:
            candidate = z
        elif candidate != z:
            return f_good, g_good, "none", None
    if all_slope:
        return f_good, g_good, "all", None
    return f_good, g_good, "one", candidate


def bad_slopes_by_support_linear(
    field: raw.Field,
    domain: list[tuple[int, ...]],
    k: int,
    agreement: int,
    f_values: dict[tuple[int, ...], tuple[int, ...]],
    g_values: dict[tuple[int, ...], tuple[int, ...]],
) -> dict[str, Any]:
    slopes: set[tuple[int, ...]] = set()
    all_slopes = False
    support_count = 0
    skipped_same_set = 0
    all_slope_supports = 0
    one_slope_supports = 0
    none_supports = 0
    for size in range(agreement, len(domain) + 1):
        for support in combinations(domain, size):
            support_count += 1
            f_good, g_good, mode, slope = support_good_and_slope(
                field, f_values, g_values, support, k
            )
            if f_good and g_good:
                skipped_same_set += 1
                continue
            if mode == "all":
                all_slopes = True
                all_slope_supports += 1
            elif mode == "one":
                one_slope_supports += 1
                assert slope is not None
                slopes.add(slope)
            else:
                none_supports += 1
    if all_slopes:
        slopes = set(field.elements())
    return {
        "bad_slopes": slopes,
        "support_stats": {
            "supports_checked": support_count,
            "skipped_same_set_supports": skipped_same_set,
            "all_slope_supports": all_slope_supports,
            "one_slope_supports": one_slope_supports,
            "none_supports": none_supports,
        },
    }


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


def first_full_degree_alpha(field: raw.Field) -> tuple[int, ...]:
    for x in field.elements():
        if x != field.zero and raw.minimal_field_degree(field, x) == field.degree:
            return x
    raise AssertionError("no full-degree alpha")


def scan_prime(p: int) -> dict[str, Any]:
    extension_degree = 4
    field = raw.Field(p, extension_degree, raw.find_irreducible(p, extension_degree))
    domain = [field.embed(x) for x in range(1, p)]
    n = p - 1
    k = n // 2
    agreement = frontier_agreement(n, k)
    alpha = first_full_degree_alpha(field)
    f0, g0 = raw.candidate_values(field, "pole_conjugate_pair", alpha, domain)
    component_sets: list[set[tuple[int, ...]]] = []
    component_reports = []
    support_stats = []
    for step in range(extension_degree):
        f_i = conjugate_values(field, f0, step)
        g_i = conjugate_values(field, g0, step)
        solved = bad_slopes_by_support_linear(field, domain, k, agreement, f_i, g_i)
        slopes = solved["bad_slopes"]
        component_sets.append(slopes)
        support_stats.append(solved["support_stats"])
        expected = {field.frobenius(slope, step) for slope in component_sets[0]}
        component_reports.append(
            {
                "component": step,
                "equals_frobenius_conjugate_of_component_0": slopes == expected,
                **slope_payload(field, slopes),
            }
        )

    union_slopes = set().union(*component_sets)
    intersection_slopes = set.intersection(*component_sets)
    union_payload = slope_payload(field, union_slopes)
    complete_orbits = union_payload["frobenius_orbit_diagnostic"]["complete_full_orbits"]
    return {
        "p": p,
        "extension_degree": extension_degree,
        "field_size": field.q,
        "irreducible_modulus_coefficients_low_to_high": list(field.modulus),
        "domain": "F_p^*",
        "n": n,
        "k": k,
        "agreement": agreement,
        "alpha": field.encode(alpha),
        "support_stats_by_component": support_stats,
        "component_reports": component_reports,
        "union": union_payload,
        "intersection": slope_payload(field, intersection_slopes),
        "complete_K_eq_F_orbits": complete_orbits,
        "complete_K_eq_F_orbits_per_p": f"{complete_orbits}/{p}",
        "status": "PASS",
    }


def build_certificate() -> dict[str, Any]:
    rows = [scan_prime(p) for p in (7, 11)]
    cert = {
        "schema": "f1-extension-quartic-growth-v1",
        "status": "EXPERIMENTAL_AUDIT",
        "source_packet": "f1-extension-galois-packaging-v1",
        "target": (
            "Test whether the quartic pole-conjugate-pair K=F orbit package "
            "stays bounded or grows with the base prime."
        ),
        "rows": rows,
        "summary": {
            "base_primes": [row["p"] for row in rows],
            "complete_K_eq_F_orbits": [
                row["complete_K_eq_F_orbits"] for row in rows
            ],
            "full_degree_union_slopes": [
                row["union"]["full_degree_bad_slopes"] for row in rows
            ],
            "growth_reading": (
            "The quartic pole-conjugate-pair package grows from 14 complete "
            "K=F orbits at p=7 to 166 complete K=F orbits at p=11.  This "
            "is not bounded at the first larger row and remains a live "
            "extension-cell obstruction candidate."
            ),
        },
        "non_claims": [
            "does not prove asymptotic growth",
            "does not classify all quartic F-valued pairs",
            "does not prove a deployed paid_extension counterexample",
            "does not close the descended-cycle inventory problem",
        ],
    }
    validate(cert)
    return cert


def validate(cert: dict[str, Any]) -> None:
    if cert.get("schema") != "f1-extension-quartic-growth-v1":
        raise AssertionError("unexpected schema")
    if cert.get("status") != "EXPERIMENTAL_AUDIT":
        raise AssertionError("unexpected status")
    expected = {
        7: {"orbits": 14, "full": 56, "total": 58, "buckets": {"1": 0, "2": 2, "4": 56}},
        11: {
            "orbits": 166,
            "full": 664,
            "total": 668,
            "buckets": {"1": 0, "2": 4, "4": 664},
        },
    }
    for row in cert["rows"]:
        e = expected[row["p"]]
        if row["complete_K_eq_F_orbits"] != e["orbits"]:
            raise AssertionError(f"unexpected orbit count for p={row['p']}")
        if row["union"]["full_degree_bad_slopes"] != e["full"]:
            raise AssertionError(f"unexpected full-degree count for p={row['p']}")
        if row["union"]["total_bad_slopes"] != e["total"]:
            raise AssertionError(f"unexpected total count for p={row['p']}")
        if row["union"]["minimal_field_degree_counts"] != e["buckets"]:
            raise AssertionError(f"unexpected buckets for p={row['p']}")
        diag = row["union"]["frobenius_orbit_diagnostic"]
        if not diag["full_orbit_closed_for_fixed_line"]:
            raise AssertionError(f"union not orbit-closed for p={row['p']}")
        if row["intersection"]["total_bad_slopes"] != 0:
            raise AssertionError(f"unexpected common-slope intersection for p={row['p']}")
        if not all(
            report["equals_frobenius_conjugate_of_component_0"]
            for report in row["component_reports"]
        ):
            raise AssertionError(f"component conjugacy failed for p={row['p']}")
    if cert["summary"]["complete_K_eq_F_orbits"] != [14, 166]:
        raise AssertionError("summary orbit counts changed")


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
    bad["rows"][1]["complete_K_eq_F_orbits"] -= 1
    try:
        validate(bad)
    except AssertionError:
        return
    raise AssertionError("tamper selftest failed")


def print_summary(cert: dict[str, Any]) -> None:
    print("f1_extension_quartic_growth")
    print(f"  status: {cert['status']}")
    print(f"  reading: {cert['summary']['growth_reading']}")
    for row in cert["rows"]:
        print(
            f"  F_{row['p']}^4/F_{row['p']}: n={row['n']} a={row['agreement']} "
            f"union={row['union']['total_bad_slopes']} "
            f"K=F={row['union']['full_degree_bad_slopes']} "
            f"orbits={row['complete_K_eq_F_orbits']}"
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
