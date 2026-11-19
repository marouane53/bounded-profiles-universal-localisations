"""Regression checks for the JSON family tool, including empty and nonsplit tuples."""
from itertools import product
from random import Random

from check_lib import require_asserts
from compare_operator_pairs import compare
from family_tools import compute


def scalars(values):
    return [[[v]] for v in values]


def scalar_orbit(w, n, q, p):
    # Compute each coefficient directly from the quotient relation, rather
    # than by repeated generator substitution. Here Q is one-dimensional.
    return {tuple(w[(i+j) % n] * pow(q, (i+j)//n, p) % p for i in range(n))
            for j in range(n*(p-1))}


def main():
    require_asserts()
    count = 0
    for p, n, q in ((2,3,1),(3,2,2),(5,2,4)):
        points = list(product(range(p), repeat=n))
        for a, b in product(points, repeat=2):
            ab = compute(dict(operation="compare", prime=p, n=n, Q=[[q]],
                              first=scalars(a), second=scalars(b)))
            same = b in scalar_orbit(a, n, q, p)
            assert ab["equal"] == ab["first_subset_of_second"] == ab["second_subset_of_first"] == same
            count += 1
    rng = Random(39)
    for p in (2,3,5):
        for _ in range(25):
            def pair(d):
                return [[[rng.randrange(p) for _ in range(d)] for _ in range(d)] for _ in range(2)]
            params = dict(prime=p, first=pair(rng.randrange(3)), second=pair(rng.randrange(3)))
            old = compare(params)
            new = compute(dict(params, operation="compare", n=2, Q=[[1]]))
            assert all(new[k] == val for k, val in old.items())
    # Different dimensions, an extension with the same simple factor, and a strict inclusion.
    zero = scalars([0,0,0,0])
    nonsplit = [[[0,1],[0,0]]] * 4
    both = [[[0,0],[0,1]]] * 4
    base = dict(operation="compare", prime=3, n=2, Q=[[1,0],[0,1]], second=zero)
    assert compute(dict(base, first=nonsplit))["equal"]
    proper = compute(dict(base, first=both))
    assert proper["first_subset_of_second"] and not proper["second_subset_of_first"]
    assert compute(dict(base, first=[[]] * 4))["second_subset_of_first"]
    c = compute(dict(operation="cubic", prime=3, C=[[1,2],[2,1]], degree=5, direct_bar=True))
    assert c["rank"] == 1 and c["direct_bar_dimensions"] == c["residue_field_ext_dimensions"]
    w = compute(dict(operation="word", prime=5, n=3, Q=[[1,2],[0,1]], word=[0,0,0,1]))
    assert w["normal_form"] == [{"word": (1,0,0,0), "coefficient": 1},
                                {"word": (2,0,0,0), "coefficient": 2}]
    print(f"Family tool matched {count} scalar comparisons and 75 comparisons against the earlier pair tool.")
    print("Checked multi-generator nonsplit and empty tuples, strict inclusion, direct bar and a mixed normal form.")


if __name__ == "__main__":
    main()
