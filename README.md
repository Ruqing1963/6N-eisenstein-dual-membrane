# The Dual Membrane and Korselt Colour-Blindness in Eisenstein Carmichael Ideals

**Part XXXVIII — Arithmetic Geodynamics on the 6N Skeleton**

Ruqing Chen · GUT Geoservice Inc., Montreal, Quebec, Canada · ruqing@hotmail.com

---

This is the closing part of the *Arithmetic Geodynamics on the 6N Skeleton* series. It lifts the
one-dimensional 6N construction — the integers with the primes 2 and 3 stripped, every surviving prime
on a wing 6N±1 — to its two-dimensional source, the Eisenstein integers **Z[ω]** (ω = e^{2πi/3}), and
studies Carmichael ideals there.

## Results

**1. The dual membrane theorem (a rigid exterior boundary).**
No Eisenstein Carmichael ideal is divisible by the prime above 3 or the prime above 2:

- **Ramified membrane (3-adic).** The ramified prime (1−ω), norm 3, divides no Eisenstein Carmichael
  ideal. Every non-ramified Eisenstein prime has norm ≡ 1 (mod 3), so each forces 3 | N(α)−1; the
  ramified prime would force N(α) ≡ 0 (mod 3) — contradiction.
- **Inert-2 membrane (2-adic).** The inert prime 2, norm 4, divides no Eisenstein Carmichael ideal.
  Its presence makes N(α) ≡ 0 (mod 4), so N(α)−1 is odd, while every other factor contributes an even
  N(π)−1 that cannot divide it.

The two exclusions rest on different valuations and are logically independent. Together they show that
the removal of 2 and 3 — imposed *by definition* in one dimension — is in two dimensions a **theorem**:
the Carmichael structure expels the primes above 2 and 3 on its own. The 6N skeleton ceases to be an
axiom and becomes a consequence of the Korselt criterion.

**2. Korselt colour-blindness (a free interior — a measured negative result).**
Above the membrane the surviving factors are of two unit-invariant types — split (norm p ≡ 1 mod 3) and
higher inert (norm q², q ≡ 2 mod 3, q ≥ 5). The type-composition spectrum of Eisenstein Carmichael
ideals is statistically indistinguishable from the null model induced by the factor pool alone:

| | pure-split | mixed | pure-inert |
|---|---|---|---|
| **Eisenstein Carmichael** | 0.663 | 0.337 | 0.000 |
| **Null (norm ceiling only, no Korselt)** | 0.642 | 0.354 | 0.004 |

Within the resolution measured, the Korselt criterion is **colour-blind**: it constrains only the
multiplicative congruence of norms and expresses no preference among factor types. The one-dimensional
monochromatic collapse has no type-level analogue in Z[ω].

## Repository contents

```
eisenstein_dual_membrane.py            reproducer (standard library only)
data/eisenstein_carmichael_examples.csv enumerated Eisenstein Carmichael ideals
                                        (factor norms + split/inert type labels)
figures/eisenstein_membrane_figures.*   (A) instrument gate occupancy; (B) colour-blindness spectrum
paper/Chen_6N_Paper38.tex / .pdf        the paper
```

## Reproducing the results

No external dependencies (Python 3, standard library only):

```bash
python3 eisenstein_dual_membrane.py
```

This regenerates, from scratch:

- **[GATE]** the instrument gate — Eisenstein prime sextant occupancy is exactly 1/6 at
  N ≤ 4×10⁴, 10⁵, 2.5×10⁵ (max deviation 0.000000), confirming the sieve is unbiased;
- **[MEMBRANE]** the dual-membrane confirmation — among the 95 Eisenstein Carmichael ideals with
  N(α) ≤ 10⁷ from a factor pool that explicitly **includes** the ramified norm 3 and the inert norm 4,
  exactly zero use them; every ideal has N(α) ≡ 1 (mod 3) and N(α) odd;
- **[COLOUR]** the colour-blindness comparison — pure-split fraction 0.663 (Carmichael) vs 0.642
  (null), agreeing to ~2%.

The figure (`figures/`) was produced with `matplotlib`; the core reproducer needs no third-party
packages.

## A note on what is and is not new

In keeping with the honest-ledger discipline of this series: Carmichael ideals over rings of integers
are a defined and studied object (Steele 2008; Bae–Hu–Sha 2023), and the dual membrane is most likely
an elementary corollary of that abstract theory rather than a new algebraic theorem. The contribution
claimed here is its *placement* — the recognition that the 6N skeleton's defining exclusion of 2 and 3
is exactly the 2-adic and 3-adic geometric blockage of the primes above 2 and 3 in the six-fold phase
space of Z[ω] — together with the first large-scale empirical map of the Carmichael type spectrum above
the membrane, which establishes the colour-blindness against an explicit null model.

## Citation

See `CITATION.cff`. The series' Part XXXVI (rational Carmichael collapse) is archived at
[doi:10.5281/zenodo.20611583](https://doi.org/10.5281/zenodo.20611583).

## License

Code under MIT (`LICENSE`). The paper and dataset are released for open scholarly use; please cite.
