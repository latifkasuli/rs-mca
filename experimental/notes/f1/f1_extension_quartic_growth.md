# F1 extension quartic growth scan

Status: EXPERIMENTAL / AUDIT.

This packet follows the `F_7^4/F_7` Galois-packaged
`pole_conjugate_pair` signal from `f1_extension_galois_packaging.md` to the
next quartic toy row.  It specializes the support-wise MCA scan to the quartic
family and solves the per-support slope constraint directly, so the
`F_11^4/F_11` row is still exact at toy scale.

## Result

Frozen certificate:
`experimental/data/certificates/f1-extension-quartic-growth/f1_extension_quartic_growth.json`.

Replay:

```bash
python3 experimental/scripts/verify_f1_extension_quartic_growth.py --check \
  experimental/data/certificates/f1-extension-quartic-growth/f1_extension_quartic_growth.json
```

For the Frobenius-packaged quartic `pole_conjugate_pair` family:

| row | `n` | `k` | `a` | union total slopes | `K=F` slopes | complete `K=F` Frobenius orbits |
|---|---:|---:|---:|---:|---:|---:|
| `F_7^4/F_7` | 6 | 3 | 4 | 58 | 56 | 14 |
| `F_11^4/F_11` | 10 | 5 | 6 | 668 | 664 | 166 |

The `F_11^4/F_11` row has minimal-field buckets
`{1:0, 2:4, 4:664}` in the packaged union.  Each individual conjugate
component has `207` slopes, with buckets `{1:0, 2:2, 4:205}`; the full
Galois package is orbit-closed and has empty common intersection across all
four conjugate components.

## Interpretation

The quartic package does **not** stay bounded at the first larger base-prime
row.  The number of complete full-degree Frobenius orbits grows from `14` at
`p=7` to `166` at `p=11`.

This keeps the family alive as an extension-cell obstruction candidate.  The
next proof-level question is not whether the first toy packet was a finite
accident; it is how the descended-cycle inventory classifies this growing
quartic package, or whether a paid-extension branch must charge it explicitly.

## Non-claims

- No asymptotic growth theorem.
- No classification of all quartic `F`-valued received lines.
- No proof that the family survives deployed extension-cell ledgers.
- No adjacent safe upper certificate.
- No proof that the observed growth is not absorbed by a descended-cycle
  classification rule.

## Next step

Optimize or derive the descended-cycle classification for this quartic package.
If moving toward proof, useful references are Weil restriction, finite Galois
descent for varieties and cycles, and Lang-Weil style finite-field point
bounds.  If staying computational, the next datum is a further specialized row
or a symbolic description of the `p=11` slope-support constraints.
