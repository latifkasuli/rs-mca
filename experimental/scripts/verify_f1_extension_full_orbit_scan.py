#!/usr/bin/env python3
"""Toy scan for the F1 extension full-orbit branch.

The frontier-adjacent extension-cell packet leaves one branch open: genuinely
F-valued bad slopes whose minimal field is the whole extension K=F.  This
script does a small exact falsifier scan on tower analogues F_{p^e}/F_p.  It is
not an extension-cell theorem and not a deployed-row certificate.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from itertools import combinations, product
from math import ceil
from pathlib import Path
from typing import Any, Iterable


REPO = Path(__file__).resolve().parents[2]
ARTIFACT = (
    REPO
    / "experimental"
    / "data"
    / "certificates"
    / "f1-extension-full-orbit-scan"
    / "f1_extension_full_orbit_scan.json"
)


def prime_factors(n: int) -> list[int]:
    factors: list[int] = []
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors.append(d)
            n //= d
        d += 1
    if n > 1:
        factors.append(n)
    return factors


def divisors(n: int) -> list[int]:
    return [d for d in range(1, n + 1) if n % d == 0]


class Field:
    def __init__(self, p: int, degree: int, modulus: tuple[int, ...]):
        self.p = p
        self.degree = degree
        self.modulus = modulus
        self.q = p**degree
        self.zero = (0,) * degree
        self.one = (1,) + (0,) * (degree - 1)

    def embed(self, x: int) -> tuple[int, ...]:
        return (x % self.p,) + (0,) * (self.degree - 1)

    def add(self, a: tuple[int, ...], b: tuple[int, ...]) -> tuple[int, ...]:
        return tuple((x + y) % self.p for x, y in zip(a, b))

    def sub(self, a: tuple[int, ...], b: tuple[int, ...]) -> tuple[int, ...]:
        return tuple((x - y) % self.p for x, y in zip(a, b))

    def neg(self, a: tuple[int, ...]) -> tuple[int, ...]:
        return tuple((-x) % self.p for x in a)

    def mul(self, a: tuple[int, ...], b: tuple[int, ...]) -> tuple[int, ...]:
        tmp = [0] * (2 * self.degree - 1)
        for i, ai in enumerate(a):
            if ai == 0:
                continue
            for j, bj in enumerate(b):
                if bj:
                    tmp[i + j] = (tmp[i + j] + ai * bj) % self.p
        for i in range(len(tmp) - 1, self.degree - 1, -1):
            coeff = tmp[i] % self.p
            if coeff == 0:
                continue
            offset = i - self.degree
            for j, mj in enumerate(self.modulus):
                tmp[offset + j] = (tmp[offset + j] - coeff * mj) % self.p
        return tuple(x % self.p for x in tmp[: self.degree])

    def pow(self, a: tuple[int, ...], exponent: int) -> tuple[int, ...]:
        if exponent < 0:
            return self.pow(self.inv(a), -exponent)
        out = self.one
        cur = a
        e = exponent
        while e:
            if e & 1:
                out = self.mul(out, cur)
            cur = self.mul(cur, cur)
            e >>= 1
        return out

    def inv(self, a: tuple[int, ...]) -> tuple[int, ...]:
        if a == self.zero:
            raise ZeroDivisionError("field inverse of zero")
        return self.pow(a, self.q - 2)

    def div(self, a: tuple[int, ...], b: tuple[int, ...]) -> tuple[int, ...]:
        return self.mul(a, self.inv(b))

    def frobenius(self, a: tuple[int, ...], steps: int = 1) -> tuple[int, ...]:
        out = a
        for _ in range(steps):
            out = self.pow(out, self.p)
        return out

    def elements(self) -> Iterable[tuple[int, ...]]:
        yield from product(range(self.p), repeat=self.degree)

    def encode(self, a: tuple[int, ...]) -> list[int]:
        return list(a)


def poly_mulmod(a: list[int], b: list[int], modulus: tuple[int, ...], p: int) -> list[int]:
    degree = len(modulus) - 1
    tmp = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            tmp[i + j] = (tmp[i + j] + ai * bj) % p
    for i in range(len(tmp) - 1, degree - 1, -1):
        coeff = tmp[i] % p
        if coeff == 0:
            continue
        offset = i - degree
        for j, mj in enumerate(modulus):
            tmp[offset + j] = (tmp[offset + j] - coeff * mj) % p
    return (tmp[:degree] + [0] * degree)[:degree]


def poly_powmod(base: list[int], exponent: int, modulus: tuple[int, ...], p: int) -> list[int]:
    degree = len(modulus) - 1
    out = [1] + [0] * (degree - 1)
    cur = base[:degree]
    e = exponent
    while e:
        if e & 1:
            out = poly_mulmod(out, cur, modulus, p)
        cur = poly_mulmod(cur, cur, modulus, p)
        e >>= 1
    return out


def is_irreducible(modulus: tuple[int, ...], p: int) -> bool:
    degree = len(modulus) - 1
    x = [0, 1] + [0] * (degree - 2)
    if poly_powmod(x, p**degree, modulus, p) != x:
        return False
    for prime in set(prime_factors(degree)):
        test = poly_powmod(x, p ** (degree // prime), modulus, p)
        if test == x:
            return False
    return True


def find_irreducible(p: int, degree: int) -> tuple[int, ...]:
    for coeffs in product(range(p), repeat=degree):
        if coeffs[0] == 0:
            continue
        modulus = tuple(coeffs) + (1,)
        if is_irreducible(modulus, p):
            return modulus
    raise RuntimeError(f"no irreducible polynomial found over F_{p} of degree {degree}")


def interpolate(field: Field, points: list[tuple[tuple[int, ...], tuple[int, ...]]]) -> list[tuple[int, ...]]:
    coeffs: list[tuple[int, ...]] = []
    for i, (xi, yi) in enumerate(points):
        basis = [field.one]
        denom = field.one
        for j, (xj, _yj) in enumerate(points):
            if i == j:
                continue
            nxt = [field.zero] * (len(basis) + 1)
            for d, bd in enumerate(basis):
                nxt[d] = field.sub(nxt[d], field.mul(xj, bd))
                nxt[d + 1] = field.add(nxt[d + 1], bd)
            basis = nxt
            denom = field.mul(denom, field.sub(xi, xj))
        scale = field.div(yi, denom)
        while len(coeffs) < len(basis):
            coeffs.append(field.zero)
        for d, bd in enumerate(basis):
            coeffs[d] = field.add(coeffs[d], field.mul(scale, bd))
    while coeffs and coeffs[-1] == field.zero:
        coeffs.pop()
    return coeffs


def poly_eval(field: Field, coeffs: list[tuple[int, ...]], x: tuple[int, ...]) -> tuple[int, ...]:
    out = field.zero
    for coeff in reversed(coeffs):
        out = field.add(field.mul(out, x), coeff)
    return out


def degree_lt_k_on_support(
    field: Field,
    domain: list[tuple[int, ...]],
    values: dict[tuple[int, ...], tuple[int, ...]],
    support: tuple[tuple[int, ...], ...],
    k: int,
) -> bool:
    if len(support) <= k:
        return True
    base = list(support[:k])
    coeffs = interpolate(field, [(x, values[x]) for x in base])
    if len(coeffs) > k:
        return False
    return all(poly_eval(field, coeffs, x) == values[x] for x in support)


def minimal_field_degree(field: Field, x: tuple[int, ...]) -> int:
    for d in divisors(field.degree):
        if field.frobenius(x, d) == x:
            return d
    raise AssertionError("Frobenius orbit did not close")


def full_orbit_diagnostic(field: Field, slopes: list[tuple[int, ...]]) -> dict[str, Any]:
    """Measure Frobenius closure of full-degree slopes for one fixed line.

    A raw F-valued received line is not itself Galois-stable in general, so the
    bad slopes of that fixed line need not be a union of Frobenius orbits.  The
    extension-cell descended-cycle target is orbit-shaped after descent; this
    diagnostic records whether the sampled fixed-line data already has that
    stronger property instead of assuming it.
    """

    slope_set = set(slopes)
    seen: set[tuple[int, ...]] = set()
    touched = []
    complete = 0
    for slope in sorted(slopes):
        if slope in seen:
            continue
        orbit = []
        cur = slope
        for _ in range(field.degree):
            orbit.append(cur)
            cur = field.frobenius(cur)
        orbit_set = set(orbit)
        if cur != slope or len(orbit_set) != field.degree:
            raise AssertionError("bad full Frobenius orbit")
        present = sorted(orbit_set & slope_set)
        seen |= orbit_set
        is_complete = len(present) == field.degree
        complete += int(is_complete)
        touched.append(
            {
                "present": [field.encode(x) for x in present],
                "missing_count": field.degree - len(present),
                "complete": is_complete,
            }
        )
    return {
        "full_orbit_closed_for_fixed_line": all(item["complete"] for item in touched),
        "touched_full_orbits": len(touched),
        "complete_full_orbits": complete,
        "orbits": touched,
    }


def candidate_values(
    field: Field,
    kind: str,
    alpha: tuple[int, ...],
    domain: list[tuple[int, ...]],
) -> tuple[dict[tuple[int, ...], tuple[int, ...]], dict[tuple[int, ...], tuple[int, ...]]]:
    f_values: dict[tuple[int, ...], tuple[int, ...]] = {}
    g_values: dict[tuple[int, ...], tuple[int, ...]] = {}
    alpha_p = field.frobenius(alpha)
    alpha2 = field.mul(alpha, alpha)
    for x in domain:
        inv = field.inv(field.sub(x, alpha))
        if kind == "pole_const_linear":
            f, g = inv, field.mul(x, inv)
        elif kind == "pole_alpha_const":
            f, g = field.mul(alpha, inv), inv
        elif kind == "pole_conjugate_pair":
            inv_p = field.inv(field.sub(x, alpha_p))
            f, g = inv, inv_p
        elif kind == "mixed_quadratic":
            x2 = field.mul(x, x)
            f = field.add(field.mul(alpha, x2), field.mul(alpha_p, x))
            g = field.add(field.mul(alpha2, x), alpha)
        else:
            raise ValueError(f"unknown candidate kind {kind}")
        f_values[x] = f
        g_values[x] = g
    return f_values, g_values


def exact_bad_slope_count(
    field: Field,
    domain: list[tuple[int, ...]],
    k: int,
    agreement: int,
    f_values: dict[tuple[int, ...], tuple[int, ...]],
    g_values: dict[tuple[int, ...], tuple[int, ...]],
) -> dict[str, Any]:
    supports = [
        tuple(support)
        for size in range(agreement, len(domain) + 1)
        for support in combinations(domain, size)
    ]
    f_good = {
        support: degree_lt_k_on_support(field, domain, f_values, support, k)
        for support in supports
    }
    g_good = {
        support: degree_lt_k_on_support(field, domain, g_values, support, k)
        for support in supports
    }
    bad_slopes: dict[tuple[int, ...], dict[str, Any]] = {}
    for slope in field.elements():
        h_values = {
            x: field.add(f_values[x], field.mul(slope, g_values[x]))
            for x in domain
        }
        for support in supports:
            if f_good[support] and g_good[support]:
                continue
            if degree_lt_k_on_support(field, domain, h_values, support, k):
                bad_slopes[slope] = {
                    "support_size": len(support),
                    "support_base_values": [x[0] for x in support],
                    "f_good": f_good[support],
                    "g_good": g_good[support],
                }
                break

    degree_buckets: dict[str, int] = {str(d): 0 for d in divisors(field.degree)}
    for slope in bad_slopes:
        degree_buckets[str(minimal_field_degree(field, slope))] += 1
    full_degree_slopes = sorted(
        [slope for slope in bad_slopes if minimal_field_degree(field, slope) == field.degree]
    )
    orbit_diagnostic = full_orbit_diagnostic(field, full_degree_slopes)
    return {
        "total_bad_slopes": len(bad_slopes),
        "minimal_field_degree_counts": degree_buckets,
        "full_orbit_bad_slopes": len(full_degree_slopes),
        "full_orbit_divisible": len(full_degree_slopes) % field.degree == 0,
        "frobenius_orbit_diagnostic": orbit_diagnostic,
        "examples": [
            {
                "slope": field.encode(slope),
                **bad_slopes[slope],
            }
            for slope in sorted(bad_slopes)[:6]
        ],
    }


@dataclass(frozen=True)
class RowConfig:
    p: int
    extension_degree: int
    alpha_limit: int = 2


def row_agreement(n: int, k: int) -> int:
    return max(k + 1, ceil((1116048 / 2097152) * n))


def scan_row(config: RowConfig) -> dict[str, Any]:
    p, degree = config.p, config.extension_degree
    if p <= 2:
        return {
            "p": p,
            "extension_degree": degree,
            "status": "SKIPPED_NO_NONTRIVIAL_MULTIPLICATIVE_DOMAIN",
            "reason": "F_p^* has order <= 1, so there is no rate-1/2 smooth multiplicative toy row",
        }
    n = p - 1
    if n % 2:
        return {
            "p": p,
            "extension_degree": degree,
            "status": "SKIPPED_ODD_DOMAIN_ORDER",
            "reason": "using D=F_p^*, rate 1/2 requires even n",
        }
    k = n // 2
    agreement = row_agreement(n, k)
    modulus = find_irreducible(p, degree)
    field = Field(p, degree, modulus)
    domain = [field.embed(x) for x in range(1, p)]
    alphas = [
        x
        for x in field.elements()
        if x != field.zero and minimal_field_degree(field, x) == degree
    ][: config.alpha_limit]
    kinds = [
        "pole_const_linear",
        "pole_alpha_const",
        "pole_conjugate_pair",
        "mixed_quadratic",
    ]

    candidates = []
    for alpha_index, alpha in enumerate(alphas):
        for kind in kinds:
            f_values, g_values = candidate_values(field, kind, alpha, domain)
            count = exact_bad_slope_count(field, domain, k, agreement, f_values, g_values)
            candidates.append(
                {
                    "case": f"p{p}_e{degree}_alpha{alpha_index}_{kind}",
                    "kind": kind,
                    "alpha": field.encode(alpha),
                    "f_is_codeword": degree_lt_k_on_support(
                        field, domain, f_values, tuple(domain), k
                    ),
                    "g_is_codeword": degree_lt_k_on_support(
                        field, domain, g_values, tuple(domain), k
                    ),
                    **count,
                }
            )
    max_full = max((row["full_orbit_bad_slopes"] for row in candidates), default=0)
    max_total = max((row["total_bad_slopes"] for row in candidates), default=0)
    return {
        "p": p,
        "extension_degree": degree,
        "field_size": field.q,
        "irreducible_modulus_coefficients_low_to_high": list(modulus),
        "domain": "F_p^*",
        "n": n,
        "k": k,
        "agreement": agreement,
        "candidate_count": len(candidates),
        "max_total_bad_slopes": max_total,
        "max_full_orbit_bad_slopes": max_full,
        "max_full_orbit_exceeds_one_orbit": max_full > degree,
        "status": "PASS",
        "candidates": candidates,
    }


def build_certificate() -> dict[str, Any]:
    row_configs = [
        RowConfig(2, 4),
        RowConfig(3, 4),
        RowConfig(5, 4),
        RowConfig(7, 4),
        RowConfig(3, 6),
        RowConfig(5, 6),
    ]
    rows = [scan_row(config) for config in row_configs]
    scanned = [row for row in rows if row["status"] == "PASS"]
    full_counts = [row["max_full_orbit_bad_slopes"] for row in scanned]
    max_full = max(full_counts, default=0)
    p_scale_signals = [
        {
            "p": row["p"],
            "extension_degree": row["extension_degree"],
            "max_full_orbit_bad_slopes": row["max_full_orbit_bad_slopes"],
            "responsible_cases": [
                candidate["case"]
                for candidate in row["candidates"]
                if candidate["full_orbit_bad_slopes"] == row["max_full_orbit_bad_slopes"]
            ],
            "all_responsible_cases_frobenius_closed": all(
                candidate["frobenius_orbit_diagnostic"]["full_orbit_closed_for_fixed_line"]
                for candidate in row["candidates"]
                if candidate["full_orbit_bad_slopes"] == row["max_full_orbit_bad_slopes"]
            ),
        }
        for row in scanned
        if row["max_full_orbit_bad_slopes"] >= row["p"]
    ]
    cert = {
        "schema": "f1-extension-full-orbit-scan-v1",
        "status": "EXPERIMENTAL_AUDIT",
        "upstream_target": (
            "frontier paid_extension K=F full-orbit branch; see "
            "frontier_extension_cell_targets_v1"
        ),
        "scan_scope": (
            "exact support-wise MCA scan for deterministic F-valued candidate "
            "families on tiny tower rows F_{p^e}/F_p with D=F_p^*"
        ),
        "rows": rows,
        "summary": {
            "scanned_rows": len(scanned),
            "skipped_rows": len(rows) - len(scanned),
            "max_full_orbit_bad_slopes": max_full,
            "p_scale_signal_rows": p_scale_signals,
            "any_raw_fixed_line_p_scale_signal": bool(p_scale_signals),
            "reading": (
                "A raw fixed-line p-scale K=F signal appears in the quartic p=7 "
                "pole-conjugate-pair family.  The responsible fixed-line bad "
                "slope set is not Frobenius-closed, so this is a "
                "falsifier-shaped extension signal, not yet a descended-cycle "
                "classification counterexample or an extension-cell theorem."
            ),
        },
        "non_claims": [
            "does not construct the descended-cycle classification inventory",
            "does not classify all F-valued pairs",
            "does not pay the frontier paid_extension cell",
            "does not prove the adjacent safe upper ledger",
        ],
    }
    validate(cert)
    return cert


def validate(cert: dict[str, Any]) -> None:
    if cert.get("schema") != "f1-extension-full-orbit-scan-v1":
        raise AssertionError("unexpected schema")
    if cert.get("status") != "EXPERIMENTAL_AUDIT":
        raise AssertionError("unexpected status")
    rows = cert.get("rows")
    if not isinstance(rows, list) or not rows:
        raise AssertionError("missing rows")
    for row in rows:
        if row["status"].startswith("SKIPPED"):
            continue
        degree = row["extension_degree"]
        for candidate in row["candidates"]:
            diag = candidate["frobenius_orbit_diagnostic"]
            if diag["complete_full_orbits"] > diag["touched_full_orbits"]:
                raise AssertionError(f"bad orbit diagnostic in {candidate['case']}")
            if diag["full_orbit_closed_for_fixed_line"] != (
                candidate["full_orbit_bad_slopes"] == degree * diag["touched_full_orbits"]
            ):
                raise AssertionError(f"orbit closure mismatch in {candidate['case']}")
            bucket_total = sum(candidate["minimal_field_degree_counts"].values())
            if bucket_total != candidate["total_bad_slopes"]:
                raise AssertionError(f"bucket total mismatch in {candidate['case']}")
            if (
                candidate["minimal_field_degree_counts"].get(str(degree), 0)
                != candidate["full_orbit_bad_slopes"]
            ):
                raise AssertionError(f"full-degree bucket mismatch in {candidate['case']}")
    if "does not classify all F-valued pairs" not in cert.get("non_claims", []):
        raise AssertionError("missing non-claim")


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
    for row in bad["rows"]:
        if row["status"] == "PASS":
            row["candidates"][0]["full_orbit_bad_slopes"] += 1
            break
    try:
        validate(bad)
    except AssertionError:
        return
    raise AssertionError("tamper selftest failed")


def print_summary(cert: dict[str, Any]) -> None:
    print("f1_extension_full_orbit_scan")
    print(f"  status: {cert['status']}")
    print(f"  reading: {cert['summary']['reading']}")
    for row in cert["rows"]:
        if row["status"] != "PASS":
            print(f"  F_{row['p']}^{row['extension_degree']}: {row['status']}")
            continue
        print(
            "  "
            f"F_{row['p']}^{row['extension_degree']} n={row['n']} "
            f"k={row['k']} a={row['agreement']} candidates={row['candidate_count']} "
            f"max_KeqF={row['max_full_orbit_bad_slopes']} "
            f"max_total={row['max_total_bad_slopes']}"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--emit-defaults", action="store_true", help="write the frozen certificate")
    parser.add_argument("--check", type=Path, help="check an existing certificate")
    parser.add_argument("--json", action="store_true", help="print the generated certificate as JSON")
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
