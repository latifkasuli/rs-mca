# F1 extension Galois-packaging scan

Status: EXPERIMENTAL / AUDIT.

This packet follows up the `F_7^4/F_7` raw fixed-line signal from
`f1_extension_full_orbit_scan.md`.  It takes the Frobenius orbit of the
responsible `pole_conjugate_pair` line and asks the extension-cell question in
the right shape: does the signal collapse into base/tower-paid strata, or does
it become an orbit-closed `K=F` inventory entry?

## Result

Frozen certificate:
`experimental/data/certificates/f1-extension-galois-packaging/f1_extension_galois_packaging.json`.

Replay:

```bash
python3 experimental/scripts/verify_f1_extension_galois_packaging.py --check \
  experimental/data/certificates/f1-extension-galois-packaging/f1_extension_galois_packaging.json
```

For both checked full-degree choices of `alpha` in the `F_7^4/F_7`
`pole_conjugate_pair` family:

| packaged object | total slopes | `K=F` slopes | complete `K=F` Frobenius orbits | common to all conjugates |
|---|---:|---:|---:|---:|
| Frobenius orbit union | 58 | **56** | **14** | 0 |

Each individual conjugate component has the same raw profile as the original
fixed line:

```text
total slopes = 15, minimal-field buckets = {1:0, 2:1, 4:14}.
```

Those individual fixed-line bad-slope sets are not Frobenius-closed.  Their
Galois union is Frobenius-closed:

```text
K=F slopes = 56 = 14 * 4 complete full-degree orbits.
```

## Interpretation

The raw signal **does not collapse** into base/tower-paid strata after Galois
packaging.  It becomes an orbit-closed `K=F` inventory candidate in the toy
quartic row.

This is still not a frontier-cell counterexample.  The scan is tiny
(`F_7^4/F_7`, `n=6`, `k=3`, `a=4`) and uses one deterministic candidate family.
It proves neither asymptotic growth nor deployed-size obstruction.  But it
does identify a concrete component the descended-cycle inventory must classify
or charge: the quartic pole-conjugate-pair package with 14 full-degree
Frobenius orbits.

## Non-claims

- No classification of all `F`-valued received lines.
- No proof that analogous packages grow with the base prime.
- No proof that this package survives the full extension-cell ledger at
  deployed parameters.
- No adjacent safe upper certificate.

## Next step

Move from one toy package to a growth test:

1. optimize the exact scan enough to test the same `pole_conjugate_pair` family
   at the next quartic toy row;
2. determine whether the number of complete `K=F` orbits stays bounded,
   grows like `p`, or is absorbed by an identifiable descended-cycle
   classification rule;
3. if it grows, promote the family as a candidate extension-cell obstruction;
   if it stays bounded, record it as a paid finite inventory class.
