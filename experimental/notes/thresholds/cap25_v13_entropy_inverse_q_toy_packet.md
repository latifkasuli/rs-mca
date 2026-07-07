# Entropy-inverse Q toy packet

**Data:** `experimental/data/certificates/frontier-adjacent/entropy_inverse_q_toy_packet_v1.json`.
**Verifier:** `experimental/scripts/verify_entropy_inverse_q_toy_packet.py`.
**Status:** EXPERIMENTAL / AUDIT.  This is not a finite frontier-adjacent
certificate and does not prove `prob:entropy-inverse-q`.

## Purpose

`experimental/grande_finale.tex` now isolates the remaining asymptotic Q
problem as `prob:entropy-inverse-q`: a large primitive logarithmic collision
moment should either be paid by an already-named cell or force a low-rank
additive-combinatorial structure, which the moment-curve Vandermonde rank then
kills.  The proof skeleton in `rem:entropy-inverse-skeleton` has six conceptual
steps: popular fibers, signed trades, low-support trade removal, entropy BSG,
Freiman/PFR structuralization, and slice-derivative transfer back to columns.

This packet is a small exact toy for the first three observable steps.  On a few
smooth prime-field rows it enumerates the full prefix-fiber distribution,
extracts the heaviest twist-primitive fiber, and forms every signed trade
between pairs of supports in that fiber.

The mathbank context is exactly the one named in `rem:standard-inverse-gap`:
Tao--Vu inverse Littlewood-Offord gives generalized-arithmetic-progression
structure from polynomial concentration (`Additive combinatorics`, Theorem 7.22,
p. 316), and BSG/Freiman/Green--Ruzsa/coset-progressions are the standard
small-doubling structural route (`Graph Theory and Additive Combinatorics`,
Chapter 7, especially pp. 253, 256, and 280).  The packet does **not** apply
those theorems as black boxes; the current repo target is entropy/logarithmic
moment scale, not polynomial concentration.

## What is computed

For each row, the verifier uses the smooth multiplicative subgroup `D` and the
power-sum prefix key

```text
(sum_{x in S} x, sum_{x in S} x^2, ..., sum_{x in S} x^w) in F_p^w.
```

For `w < p`, this is equivalent to the first elementary coefficient prefix by
Newton identities.  A prefix value is called twist-primitive when

```text
gcd(n, {j : key_j != 0}) = 1.
```

The verifier records:

- the complete fiber histogram;
- raw and twist-primitive maximum fibers and normalized ratio
  `R = p^w * max_count / binom(n,m)`;
- the heaviest twist-primitive fiber;
- all unordered pairs inside that fiber, written as signed trades
  `A\B - B\A`;
- exact verification that the first `w` trade moments vanish;
- low-support and intersection histograms for those trades;
- a small Vandermonde-rank proxy on the union of heavy-fiber support points
  (not on the full moment-curve columns `v_y` from
  `prop:vandermonde-kills-low-rank`).

## Headline rows

| row | max twist-primitive fiber | normalized `R` | trade pairs | min trade support | union size |
|---|---:|---:|---:|---:|---:|
| `F17x_16_m8_w3` | 5 | 1.908702 | 10 | 8 | 16 |
| `mu20_F41_m10_w2` | 133 | 1.210099 | 8,778 | 6 | 20 |
| `mu24_F97_m12_w2` | 333 | 1.158660 | 55,278 | 6 | 24 |

All listed trade pairs satisfy the first `w` moment equations exactly.  The
dense rows have many exact trades, but the simple Vandermonde-rank proxy on the
union of heavy-fiber support points remains full at widths `w` and `w+1`.  Thus
the toy sees the skeleton's popular-fiber and signed-trade objects, but not a
low-rank inverse-theorem output.

## Reading

This is a small exact packet for the inverse-theorem workbench.  Its useful
message is not "Q is proved" or "Q fails"; it is:

1. twist-primitive heavy fibers are present at toy scale;
2. their support pairs produce the exact signed trades the skeleton predicts;
3. low-support trades and low Vandermonde rank are **not** automatically forced
   in these rows;
4. the next theorem still has to explain how logarithmic-moment excess becomes
   entropy-small-doubling structure, not just pairwise moment cancellation.

## Replay

```sh
python3 experimental/scripts/verify_entropy_inverse_q_toy_packet.py --check
python3 experimental/scripts/verify_entropy_inverse_q_toy_packet.py --tamper-selftest
```

The default verifier recomputes all rows exactly and compares the frozen JSON.
It is deterministic and has no random seed.

## Non-claims

- No deployed row is instantiated.
- No paid-cell first-match removal is performed beyond the twist-primitive
  classifier.
- No upper ledger, unsafe certificate, safe certificate, or adjacent pin is
  asserted.
- No polynomial-scale inverse Littlewood-Offord theorem is used as a substitute
  for the entropy-scale inverse theorem.
