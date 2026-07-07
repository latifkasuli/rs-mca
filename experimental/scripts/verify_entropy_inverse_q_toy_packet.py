#!/usr/bin/env python3
"""Exact toy packet for prob:entropy-inverse-q.

This verifier emits/checks

    experimental/data/certificates/frontier-adjacent/
        entropy_inverse_q_toy_packet_v1.json

and the companion note

    experimental/notes/thresholds/cap25_v13_entropy_inverse_q_toy_packet.md.

Scope.  This is a finite, toy-scale diagnostic for the entropy-scale inverse-Q
program in experimental/grande_finale.tex, not a deployed certificate.  For a
few small smooth-domain prefix maps, it enumerates the full m-subset
distribution, extracts the heaviest twist-primitive fiber, builds the signed
trade population from pairs of supports in that fiber, and records whether the
inverse-theorem skeleton's predicted objects are visible:

  * popular twist-primitive fibers;
  * signed trades with the first w moments vanishing;
  * low-support trade obstructions, if present;
  * a simple Vandermonde-rank proxy for the heavy fiber's union of locators.

The packet intentionally does NOT prove prob:entropy-inverse-q, does NOT remove
first-match paid cells, and does NOT move any frontier-adjacent row.  It is a
small exact "acid test" for the skeleton in rem:entropy-inverse-skeleton.

All arithmetic is over prime fields, with domains equal to smooth
multiplicative subgroups.  Prefix keys are power sums
(sum x, sum x^2, ..., sum x^w), which are equivalent to the first elementary
coefficient prefixes in these rows because w < p.
"""
from __future__ import annotations

import argparse
import contextlib
import copy
import hashlib
import io
import json
import math
import os
import sys
import tempfile
from collections import Counter, defaultdict
from itertools import combinations
from typing import Iterable


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(SCRIPT_DIR))
JSON_PATH = os.path.join(
    REPO_ROOT,
    "experimental",
    "data",
    "certificates",
    "frontier-adjacent",
    "entropy_inverse_q_toy_packet_v1.json",
)
NOTE_PATH = os.path.join(
    REPO_ROOT,
    "experimental",
    "notes",
    "thresholds",
    "cap25_v13_entropy_inverse_q_toy_packet.md",
)

WALL_ID = "CAP25-V13-ENTROPY-INVERSE-Q-TOY"

DEFAULT_ROWS = [
    {
        "label": "F17x_16_m8_w3",
        "p": 17,
        "n": 16,
        "m": 8,
        "w": 3,
        "reason": "Poisson-boundary twist-primitive prefix row; small enough to inspect all trades exactly.",
    },
    {
        "label": "mu20_F41_m10_w2",
        "p": 41,
        "n": 20,
        "m": 10,
        "w": 2,
        "reason": "Dense-bulk prefix row with nontrivial twist-primitive heavy fibers.",
    },
    {
        "label": "mu24_F97_m12_w2",
        "p": 97,
        "n": 24,
        "m": 12,
        "w": 2,
        "reason": "Largest default dense-bulk row already used by the Gamma_r packet.",
    },
]


def primitive_root(p: int) -> int:
    phi = p - 1
    factors = []
    x = phi
    d = 2
    while d * d <= x:
        if x % d == 0:
            factors.append(d)
            while x % d == 0:
                x //= d
        d += 1
    if x > 1:
        factors.append(x)
    for g in range(2, p):
        if all(pow(g, phi // q, p) != 1 for q in factors):
            return g
    raise RuntimeError(f"no primitive root for p={p}")


def subgroup_domain(p: int, n: int) -> list[int]:
    if (p - 1) % n != 0:
        raise ValueError(f"n={n} must divide p-1={p - 1}")
    g = primitive_root(p)
    h = pow(g, (p - 1) // n, p)
    vals = []
    cur = 1
    for _ in range(n):
        vals.append(cur)
        cur = (cur * h) % p
    if len(set(vals)) != n:
        raise AssertionError("subgroup generator did not have requested order")
    return sorted(vals)


def prefix_key(support: Iterable[int], p: int, w: int) -> tuple[int, ...]:
    out = [0] * w
    for x in support:
        cur = x % p
        for j in range(w):
            if j == 0:
                power = cur
            else:
                power = (power * cur) % p
            out[j] = (out[j] + power) % p
    return tuple(out)


def twist_stabilizer_size(key: tuple[int, ...], n: int) -> int:
    active = [j + 1 for j, value in enumerate(key) if value != 0]
    if not active:
        return n
    g = n
    for j in active:
        g = math.gcd(g, j)
    return g


def is_primitive_key(key: tuple[int, ...], n: int) -> bool:
    return twist_stabilizer_size(key, n) == 1


def rank_mod_p(rows: list[list[int]], p: int) -> int:
    mat = [[v % p for v in row] for row in rows if any(v % p for v in row)]
    if not mat:
        return 0
    rank = 0
    n_rows = len(mat)
    n_cols = len(mat[0])
    for col in range(n_cols):
        pivot = None
        for r in range(rank, n_rows):
            if mat[r][col] % p:
                pivot = r
                break
        if pivot is None:
            continue
        mat[rank], mat[pivot] = mat[pivot], mat[rank]
        inv = pow(mat[rank][col], -1, p)
        mat[rank] = [(v * inv) % p for v in mat[rank]]
        for r in range(n_rows):
            if r != rank and mat[r][col] % p:
                factor = mat[r][col] % p
                mat[r] = [(mat[r][c] - factor * mat[rank][c]) % p for c in range(n_cols)]
        rank += 1
        if rank == n_rows:
            break
    return rank


def vandermonde_rank_proxy(points: list[int], p: int, width: int) -> int:
    if not points:
        return 0
    rows = []
    for j in range(width):
        rows.append([pow(x, j, p) for x in points])
    return rank_mod_p(rows, p)


def sha256_json(obj) -> str:
    payload = json.dumps(obj, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def support_to_indices(support: Iterable[int], domain: list[int]) -> list[int]:
    pos = {x: i for i, x in enumerate(domain)}
    return sorted(pos[x] for x in support)


def exact_trade_population(
    supports: list[tuple[int, ...]],
    p: int,
    w: int,
    domain: list[int],
) -> dict:
    pair_count = math.comb(len(supports), 2)
    half_weight_hist = Counter()
    total_support_hist = Counter()
    intersection_hist = Counter()
    trade_vectors = Counter()
    first_bad = None
    min_half_weight = None
    max_half_weight = 0
    min_total_support = None
    union_points = set()
    example_trades = []

    for support in supports:
        union_points.update(support)

    pos = {x: i for i, x in enumerate(domain)}
    for i, supp_a in enumerate(supports):
        set_a = set(supp_a)
        for supp_b in supports[i + 1 :]:
            set_b = set(supp_b)
            plus = sorted(set_a - set_b)
            minus = sorted(set_b - set_a)
            if len(plus) != len(minus):
                first_bad = {
                    "reason": "unbalanced_trade",
                    "support_a": support_to_indices(supp_a, domain),
                    "support_b": support_to_indices(supp_b, domain),
                }
                break
            moments = []
            for j in range(1, w + 1):
                moments.append((sum(pow(x, j, p) for x in plus) - sum(pow(x, j, p) for x in minus)) % p)
            if any(moments):
                first_bad = {
                    "reason": "nonzero_trade_moment",
                    "moments": moments,
                    "support_a": support_to_indices(supp_a, domain),
                    "support_b": support_to_indices(supp_b, domain),
                }
                break
            half = len(plus)
            total = len(plus) + len(minus)
            half_weight_hist[str(half)] += 1
            total_support_hist[str(total)] += 1
            intersection_hist[str(len(set_a & set_b))] += 1
            min_half_weight = half if min_half_weight is None else min(min_half_weight, half)
            max_half_weight = max(max_half_weight, half)
            min_total_support = total if min_total_support is None else min(min_total_support, total)
            vec = [0] * len(domain)
            for x in plus:
                vec[pos[x]] = 1
            for x in minus:
                vec[pos[x]] = -1
            trade_vectors[tuple(vec)] += 1
            if len(example_trades) < 3:
                example_trades.append(
                    {
                        "plus_indices": sorted(pos[x] for x in plus),
                        "minus_indices": sorted(pos[x] for x in minus),
                        "half_weight": half,
                    }
                )
        if first_bad is not None:
            break

    if first_bad is not None:
        raise AssertionError(f"trade verification failed: {first_bad}")

    duplicate_trade_multiplicity_max = max(trade_vectors.values()) if trade_vectors else 0
    return {
        "pair_count": pair_count,
        "all_first_w_moments_vanish": True,
        "half_weight_histogram": dict(sorted(half_weight_hist.items(), key=lambda kv: int(kv[0]))),
        "total_support_histogram": dict(sorted(total_support_hist.items(), key=lambda kv: int(kv[0]))),
        "intersection_size_histogram": dict(sorted(intersection_hist.items(), key=lambda kv: int(kv[0]))),
        "min_half_weight": min_half_weight,
        "max_half_weight": max_half_weight,
        "min_total_support": min_total_support,
        "distinct_trade_vectors": len(trade_vectors),
        "duplicate_trade_multiplicity_max": duplicate_trade_multiplicity_max,
        "union_support_size": len(union_points),
        "union_vandermonde_rank_width_w": vandermonde_rank_proxy(sorted(union_points), p, w),
        "union_vandermonde_rank_width_w_plus_1": vandermonde_rank_proxy(sorted(union_points), p, w + 1),
        "example_trades": example_trades,
    }


def summarize_histogram(hist: dict[tuple[int, ...], int], p: int, n: int, m: int, w: int) -> dict:
    total = math.comb(n, m)
    expected_cells = p**w
    mean = total / expected_cells
    primitive_items = [(key, count) for key, count in hist.items() if is_primitive_key(key, n)]
    raw_key, raw_count = max(hist.items(), key=lambda kv: (kv[1], tuple(-x for x in kv[0])))
    prim_key, prim_count = max(primitive_items, key=lambda kv: (kv[1], tuple(-x for x in kv[0])))
    count_hist = Counter(hist.values())
    prim_count_hist = Counter(count for _, count in primitive_items)
    gamma2_raw_num = sum(c * c for c in hist.values()) * expected_cells
    gamma2_den = total * total
    gamma2_prim_num = sum(c * c for _, c in primitive_items) * expected_cells
    return {
        "total_supports": total,
        "prefix_cell_count": expected_cells,
        "cells_hit": len(hist),
        "primitive_cells_hit": len(primitive_items),
        "mean_fiber_size": mean,
        "raw_max": {
            "key": list(raw_key),
            "count": raw_count,
            "normalized_R": raw_count / mean,
            "twist_stabilizer_size": twist_stabilizer_size(raw_key, n),
        },
        "primitive_max": {
            "key": list(prim_key),
            "count": prim_count,
            "normalized_R": prim_count / mean,
            "twist_stabilizer_size": twist_stabilizer_size(prim_key, n),
        },
        "count_histogram_sha256": sha256_json(sorted(count_hist.items())),
        "primitive_count_histogram_sha256": sha256_json(sorted(prim_count_hist.items())),
        "gamma2_raw": gamma2_raw_num / gamma2_den,
        "gamma2_primitive": gamma2_prim_num / gamma2_den,
    }


def compute_row(row: dict) -> dict:
    p = int(row["p"])
    n = int(row["n"])
    m = int(row["m"])
    w = int(row["w"])
    domain = subgroup_domain(p, n)
    hist: dict[tuple[int, ...], int] = defaultdict(int)
    supports_by_key: dict[tuple[int, ...], list[tuple[int, ...]]] = defaultdict(list)

    for indices in combinations(range(n), m):
        support = tuple(domain[i] for i in indices)
        key = prefix_key(support, p, w)
        hist[key] += 1
        supports_by_key[key].append(support)

    summary = summarize_histogram(hist, p, n, m, w)
    heavy_key = tuple(summary["primitive_max"]["key"])
    heavy_supports = sorted(supports_by_key[heavy_key])
    trade_population = exact_trade_population(heavy_supports, p, w, domain)
    examples = [support_to_indices(s, domain) for s in heavy_supports[:5]]
    key_moments_verified = all(prefix_key(s, p, w) == heavy_key for s in heavy_supports)
    if not key_moments_verified:
        raise AssertionError("heavy-fiber supports do not share the stored key")

    low_support_trade_present = (
        trade_population["min_total_support"] is not None
        and trade_population["min_total_support"] <= 2 * w
    )

    return {
        "label": row["label"],
        "parameters": {
            "p": p,
            "n": n,
            "m": m,
            "w": w,
            "domain": "multiplicative_subgroup",
            "domain_sha256": sha256_json(domain),
        },
        "reason": row["reason"],
        "full_fiber_distribution": summary,
        "heavy_primitive_fiber": {
            "key": list(heavy_key),
            "support_count": len(heavy_supports),
            "key_moments_verified": key_moments_verified,
            "support_examples_as_domain_indices": examples,
            "all_supports_sha256": sha256_json([support_to_indices(s, domain) for s in heavy_supports]),
        },
        "signed_trade_population": trade_population,
        "inverse_skeleton_diagnostics": {
            "popular_primitive_fiber_observed": summary["primitive_max"]["normalized_R"] > 1.0,
            "low_support_trade_present_at_width_w": low_support_trade_present,
            "heavy_union_rank_equals_width_w": trade_population["union_vandermonde_rank_width_w"] == w,
            "heavy_union_rank_equals_width_w_plus_1": trade_population["union_vandermonde_rank_width_w_plus_1"]
            == min(w + 1, trade_population["union_support_size"]),
            "interpretation": (
                "The heavy twist-primitive fiber gives exact signed trades, but the "
                "toy Vandermonde-rank proxy stays full at the measured width; "
                "this is evidence for the skeleton objects, not a proof of the "
                "entropy inverse theorem."
            ),
        },
    }


def build_packet() -> dict:
    rows = [compute_row(row) for row in DEFAULT_ROWS]
    return {
        "schema": "entropy-inverse-q-toy-packet/v1",
        "wall_id": WALL_ID,
        "status": "EXPERIMENTAL / AUDIT",
        "source_problem": "experimental/grande_finale.tex prob:entropy-inverse-q and rem:entropy-inverse-skeleton",
        "companion_note": "experimental/notes/thresholds/cap25_v13_entropy_inverse_q_toy_packet.md",
        "companion_verifier": "experimental/scripts/verify_entropy_inverse_q_toy_packet.py",
        "method": {
            "prefix_model": "power-sum prefix keys (sum x^1, ..., sum x^w) over small prime-field smooth subgroups",
            "primitive_filter": "twist stabilizer gcd(n,{j: key_j != 0}) equals 1",
            "trade_model": "all unordered pairs inside the heaviest twist-primitive fiber; plus=A\\B, minus=B\\A",
            "rank_proxy": "Vandermonde rank of the union of heavy-fiber support locators at widths w and w+1",
        },
        "non_claims": [
            "Does not prove prob:entropy-inverse-q.",
            "Does not remove every first-match paid cell; primitive means twist-primitive only.",
            "Does not instantiate deployed frontier-adjacent rows.",
            "Does not move any MCA or list threshold.",
            "Does not apply polynomial-scale inverse Littlewood-Offord as a black box.",
        ],
        "mathbank_context": [
            "Tao-Vu inverse Littlewood-Offord and Green-Ruzsa/Freiman/PFR results are templates for the skeleton.",
            "The repo's current problem requires entropy/logarithmic-moment scale rather than polynomial concentration.",
        ],
        "rows": rows,
    }


def write_json(path: str, packet: dict) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(packet, fh, indent=2, sort_keys=True)
        fh.write("\n")


def load_json(path: str) -> dict:
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def compare_packet(expected: dict, actual: dict) -> list[str]:
    if expected == actual:
        return []
    diffs = []
    expected_rows = {row["label"]: row for row in expected.get("rows", [])}
    actual_rows = {row["label"]: row for row in actual.get("rows", [])}
    if set(expected_rows) != set(actual_rows):
        diffs.append(f"row labels differ: expected {sorted(expected_rows)}, actual {sorted(actual_rows)}")
    for label in sorted(set(expected_rows) & set(actual_rows)):
        e = expected_rows[label]
        a = actual_rows[label]
        for path in [
            ("full_fiber_distribution", "raw_max", "count"),
            ("full_fiber_distribution", "primitive_max", "count"),
            ("full_fiber_distribution", "primitive_max", "normalized_R"),
            ("signed_trade_population", "pair_count"),
            ("signed_trade_population", "min_total_support"),
            ("signed_trade_population", "union_support_size"),
        ]:
            ev = e
            av = a
            for key in path:
                ev = ev[key]
                av = av[key]
            if ev != av:
                diffs.append(f"{label} {'.'.join(path)} expected {ev!r} actual {av!r}")
    if not diffs:
        diffs.append("packet differs outside headline fields")
    return diffs


def check_packet_file(expected: dict, path: str) -> list[str]:
    return compare_packet(expected, load_json(path))


def tamper_selftest(packet: dict) -> None:
    checks = []

    mutated = copy.deepcopy(packet)
    mutated["rows"][0]["signed_trade_population"]["all_first_w_moments_vanish"] = False
    checks.append(("trade moment flag", mutated))

    mutated = copy.deepcopy(packet)
    mutated["rows"][0]["full_fiber_distribution"]["primitive_max"]["count"] += 1
    checks.append(("primitive max count", mutated))

    mutated = copy.deepcopy(packet)
    mutated["non_claims"] = []
    checks.append(("non-claims", mutated))

    with tempfile.TemporaryDirectory() as tmpdir:
        for name, mutated_packet in checks:
            tampered_path = os.path.join(tmpdir, f"{name.replace(' ', '_')}.json")
            write_json(tampered_path, mutated_packet)
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                rc = main(["--json-out", tampered_path, "--check"])
            if rc == 0:
                raise AssertionError(f"tamper self-test did not catch {name}")
    print(f"tamper self-test: caught {len(checks)} mutations")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json-out", default=JSON_PATH)
    parser.add_argument("--emit-defaults", action="store_true", help="write the frozen JSON packet")
    parser.add_argument("--check", action="store_true", help="recompute and compare against --json-out")
    parser.add_argument("--tamper-selftest", action="store_true")
    args = parser.parse_args(argv)

    packet = build_packet()
    if args.tamper_selftest:
        tamper_selftest(packet)
        return 0
    if args.emit_defaults:
        write_json(args.json_out, packet)
        print(f"wrote {os.path.relpath(args.json_out, REPO_ROOT)}")
    if args.check:
        diffs = check_packet_file(packet, args.json_out)
        if diffs:
            print("CHECK FAILED", file=sys.stderr)
            for diff in diffs[:20]:
                print(f"  - {diff}", file=sys.stderr)
            return 1
        print(f"CHECK OK: {os.path.relpath(args.json_out, REPO_ROOT)}")
    if not args.emit_defaults and not args.check:
        print(json.dumps(packet, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
