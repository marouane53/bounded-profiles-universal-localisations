# Two-term dg kernels and bounded profiles of universal localisations

**Marouane Lamharzi Alaoui** · marouane@masonai.fr · Preprint, 2026

This repository contains the paper *Two-term dg kernels and bounded profiles of universal localisations* and a series of research notes, dated 18 November 2026, that continue it. The material has not been peer reviewed.

## The paper

[`Lamharzi_Alaoui_2026_bounded_profiles_universal_localisations.pdf`](Lamharzi_Alaoui_2026_bounded_profiles_universal_localisations.pdf)

Every proper dg algebra over a field with cohomology in degrees −1 and 0 is realised as the dg endomorphism algebra of a perfect kernel generator for a stably flat universal localisation of finite-dimensional algebras at one map. The target is projective on both sides, and the marked kernel embedding is standard t-exact with full finite-dimensional profile. Conversely, a standard t-exact marked realisation of a connective dg algebra by an ordinary homological ring epimorphism forces a two-term projective generator. The paper also studies the bounded profiles of these embeddings and asks, as Question 3.21, whether a higher multiplication can change the profile poset at fixed graded cohomology.

The LaTeX source is in [`source/`](source). The computational supplement, with exact-arithmetic certificates for the finite examples, is in [`source/anc/supplement/`](source/anc/supplement); see its README for how to run it.

## The research notes

The notes answer Question 3.21, classify the full profile poset of a six-dimensional example, and extend the classification to an infinite family. Each note has its PDF, its LaTeX source and a script with finite exact-arithmetic checks.

| Note | Main result |
|---|---|
| [When a cubic higher product changes a bounded profile](notes/profile_classification.pdf) | For Example 3.11 both posets have two elements, but the nonzero defect categories are not equivalent. Connected quiver examples where the profile counts differ. |
| [A cubic higher multiplication creates infinitely many bounded profiles](notes/local_profile_spectrum.pdf) | A six-dimensional commutative local dg algebra whose formal model has two profiles, while the algebra itself has an infinite antichain, over every field. Six is the smallest possible dimension. |
| [Operator pairs and all bounded profiles of a six-dimensional dg algebra](notes/operator_pair_profiles.pdf) | The complete profile poset of that algebra. Two operator pairs give the same profile exactly when they have the same simple composition factors up to exchange; the remaining profiles are invertibility domains of polynomial matrices. |
| [Bounded profiles and finite-dimensional points of universal localisations](notes/support_theory/normal_localisations.pdf) | Classification of all profiles for Koszul duals with a regular normal element, including the family x^n z_i = Σ Q_ij z_j x^n for every finite-order Q, over any field. |
| [A lower bound for derived tensor products from a double derivation](notes/support_theory/two_cuts.pdf) | The closed-fibre estimate used in the classification. |
| [Central reduction for bounded profiles with several parameters](notes/support_theory/central_reduction.pdf) | A local-to-global membership test for several central parameters. |
| [A rank criterion for a cubic higher product to change bounded profiles](notes/support_theory/cubic_rank.pdf) | In a matrix family with constant cohomology, the poset changes exactly when the matrix is square and invertible. |
| [Orbit closures on representation spaces detect bounded profiles](notes/support_theory/orbit_closures.pdf) | A geometric criterion in characteristic zero and its failure in positive characteristic. |
| [Failure of constructibility for bounded Ext at fixed graded cohomology](notes/support_theory/parameter_jump.pdf) | In a flat family, the profile poset detects whether a parameter is a root of unity. |
| [Linear-recurrence zero sets in the Ext of finite modules](notes/support_theory/recurrence_ext.pdf) | Every linear-recurrence zero set occurs, up to an initial segment, in Ext between finite modules over a radical-cube-zero algebra. |

## Running the checks

The note scripts need Python 3.9 or later and only the standard library. Run each one from its own folder, for example:

```sh
cd notes && python3 verify_operator_profiles.py
cd notes/support_theory && python3 verify_normal_localisations.py
```

The scripts in `notes/support_theory/` import helpers from `notes/`, so keep the folder layout. The scripts check finite ranges only; the general statements are proved in the notes.
