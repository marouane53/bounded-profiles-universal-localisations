# Computational supplement

For *Two-term dg kernels and bounded profiles of universal localisations*, by Marouane Lamharzi Alaoui.

These exact-arithmetic certificates check the finite examples and calculations identified in Appendix E. The general theorems are proved in the manuscript; finite computations do not establish their universal claims.

## Reproduce every certificate

Use Python 3.10 or later and SymPy 1.14.0. From this directory, either run:

```sh
uv run --with sympy==1.14.0 python run_all_checks.py
```

or install the dependency in a virtual environment:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python run_all_checks.py
```

The runner executes all fifteen distinct certificates in a temporary copy. It compares each captured output byte for byte and checks both regenerated Schur multiplication tables. It does not modify the distributed scripts, outputs, or tables. A failed assertion, a changed output, or a changed table causes a nonzero exit status. A run normally takes less than a few minutes.

The complete suite was reproduced with Python 3.12.9 and SymPy 1.14.0. The older captured outputs were originally generated with Python 3.14 for the standard-library scripts and Python 3.12 for the SymPy scripts; the current reproduction matches those outputs exactly.

## Layout and scope

- `part1_comparison/`: six certificates for the comparison, intrinsic-model, and Nakayama sections. The Nakayama Hom table and the A63 stable-kernel computation use SymPy; the characteristic-two Hom table, finite-projective duality, AHKLY side conventions, and factorisation toys use the standard library. Captured files have names ending in `_output*.txt`.
- `part2_fixed_kernel/`: five certificates in `scripts/`, with matching files in `outputs/`. They cover the dg boundary, Nakayama profiles, the seven-dimensional example, formal realisation, and the group-block obstruction. Its original `run_checks.py` remains available as a runner for this group alone.
- `part3_two_term/`: `verify_nonformal.py` and the helpers `exact.py`, `m3.py`, with captured output. This checks the nineteen-dimensional non-formal realisation in Example 3.11, including multiplication, module decompositions, the strict dg map, and the transferred higher product.
- `part4_optimality/`: two new certificates. The connective-obstruction script checks the negative-degree Frobenius pairings, Schur examples, and related finite controls. The marking-boundary script checks Example 3.20. See its README for the precise scope and the two JSON multiplication tables.

- `part5_revision/`: the coordinate-swap certificate for Example 7.8, checking the algebra involution, the exchanged module actions, and the repeating resolution step.

`THEOREM_CORRESPONDENCE.csv` maps old labels to the revised numbering. Legacy theorem numbers in unchanged captured outputs refer to the previous manuscript.

Helper copies in different groups are not counted as independent certificates. All arithmetic uses rational numbers or finite prime fields. Checks over finitely many characteristics supplement, but do not replace, the arbitrary-field proofs.

## Conventions

The comparison scripts use right modules over cyclic Nakayama algebras and cohomological cones `[P_(i+1) -> P_i]` in degrees -1, 0. Their mapping-complex differentials agree with the manuscript. The scripts for the realisation and Schur examples use the left-module conventions stated in those sections. In `check_7dim.py`, the variable `N` denotes the manuscript's test object `X = C/(y)`.

Some certificates involving infinite resolutions inspect a stated finite range or period; their general extension is proved in the manuscript. The finite Hom complexes in the optimality group are checked in their entire degree range. No homological-ring-epimorphism theorem is inferred merely from a finite list of Tor vanishings.

## Integrity and licence

`SHA256SUMS` covers every distributed file in this directory except itself. Before creating a virtual environment here, check it with:

```sh
shasum -a 256 -c SHA256SUMS
```

The MIT licence in `LICENSE` applies to these scripts, outputs, and supporting computational files. It does not license the manuscript text or TeX sources.
