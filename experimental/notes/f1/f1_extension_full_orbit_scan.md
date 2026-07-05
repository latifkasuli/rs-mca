# F1 extension full-orbit toy scan

Status: EXPERIMENTAL / AUDIT.

This packet is a bounded falsifier scan for the `K=F` full-orbit branch named
by the frontier-adjacent `paid_extension` cell.  It does not pay the extension
cell, does not construct the descended-cycle inventory, and does not prove a
safe upper ledger.

## Target

The frontier-extension-cell target reduces the genuinely extension-valued
residual to the following named input: classify the `K=F` full-orbit branch, or
find a positive-dimensional extension-valued bad-slope locus.  The latter would
be a new obstruction floor for the adjacent frontier packets, because a
positive-dimensional chart contributes at least one base-field factor.

This scan tests the falsifier direction on tiny tower analogues

```text
F_{p^e} / F_p,       e in {4,6},       D = F_p^*,       rho = 1/2.
```

For each row it evaluates deterministic `F`-valued candidate pairs, computes
all support-wise MCA-bad slopes exactly, classifies each bad slope by minimal
field degree over `F_p`, and records whether the full-degree bad slopes of the
fixed line are Frobenius-closed.

## Result

Frozen certificate:
`experimental/data/certificates/f1-extension-full-orbit-scan/f1_extension_full_orbit_scan.json`.

Replay:

```bash
python3 experimental/scripts/verify_f1_extension_full_orbit_scan.py --check \
  experimental/data/certificates/f1-extension-full-orbit-scan/f1_extension_full_orbit_scan.json
```

Summary of the default rows:

| row | `n` | `k` | `a` | candidates | max total bad slopes | max `K=F` bad slopes |
|---|---:|---:|---:|---:|---:|---:|
| `F_2^4/F_2` | - | - | - | - | - | skipped: no nontrivial `F_p^*` domain |
| `F_3^4/F_3` | 2 | 1 | 2 | 8 | 1 | 1 |
| `F_5^4/F_5` | 4 | 2 | 3 | 8 | 4 | 4 |
| `F_7^4/F_7` | 6 | 3 | 4 | 8 | 15 | **14** |
| `F_3^6/F_3` | 2 | 1 | 2 | 8 | 1 | 1 |
| `F_5^6/F_5` | 4 | 2 | 3 | 8 | 4 | 4 |

The notable row is `F_7^4/F_7`: the `pole_conjugate_pair` candidates have
`14` full-degree bad slopes.  This is a `p`-scale signal (`14 >= p`) in the
raw fixed-line scan.

## Interpretation

This is a **falsifier-shaped extension signal**, not yet a frontier-cell
counterexample.

The distinction is important.  The `K=F` branch in the extension-cell target is
phrased after Galois descent into full-orbit cycle data.  The raw fixed
`F`-valued line scanned here is not generally Galois-stable, and in the
responsible `F_7^4` candidates the full-degree bad-slope set is **not**
Frobenius-closed.  So the packet does not by itself produce the descended-cycle
inventory or refute the adjacent safe ledger.

What it does provide is a concrete next target: explain whether the
`pole_conjugate_pair` signal is absorbed by the descended-cycle classification
or survives as an orbit-closed positive-dimensional branch after the correct
Galois packaging.

## Non-claims

- No classification of all `F`-valued pairs.
- No proof that the observed signal grows asymptotically.
- No proof that the observed signal is a `paid_extension` counterexample.
- No adjacent upper certificate.

## Next step

Convert the responsible `pole_conjugate_pair` row into a Galois-packaged scan:
take the Frobenius orbit of the fixed line, form the descended cycle data
requested by the extension-cell target, and decide whether the 14 raw
full-degree slopes collapse into an already-paid base/tower branch or produce
an orbit-closed `K=F` inventory entry.
