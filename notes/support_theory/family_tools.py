"""Exact prime-field computations. Read one JSON object from a file, or from -."""
import argparse
from collections import deque
import json
import sys

from check_lib import rank
from compare_operator_pairs import Basis, Field, ideal_kills_by_iteration
from verify_cubic_rank import predicted, reduced_ranks
from verify_normal_localisations import normal_form_engine


def matrix(value, field, square=False):
    if not isinstance(value, list) or not value or not all(isinstance(row, list) for row in value):
        raise ValueError("a nonempty rectangular matrix is required")
    width = len(value[0])
    if not width or any(len(row) != width for row in value) or (square and len(value) != width):
        raise ValueError("incorrect matrix shape")
    return tuple(tuple(field.parse(x) for x in row) for row in value)


def natural(value, minimum=0):
    if type(value) != int or value < minimum:
        raise ValueError(f"expected an integer at least {minimum}")
    return value


def read_tuple(value, length, field):
    if not isinstance(value, list) or len(value) != length:
        raise ValueError(f"each chart tuple needs {length} matrices, in generator/slot order")
    if all(m == [] for m in value):
        return ((),) * length
    matrices = tuple(matrix(m, field, square=True) for m in value)
    if len({len(m) for m in matrices}) != 1:
        raise ValueError("each tuple must act on one fixed vector space")
    return matrices


def twist_orbit(w, n, q, p):
    result = []
    current = w
    size = len(w[0])
    while True:
        result.append(current)
        following = []
        for i, row in enumerate(q):
            following.extend(current[i*n+1:(i+1)*n])
            following.append(tuple(
                tuple(sum(c*current[j*n][a][b] for j, c in enumerate(row)) % p
                      for b in range(size)) for a in range(size)))
        current = tuple(following)
        if current == w:
            return result


def joint_image(tuples, field):
    sizes = [len(w[0]) for w in tuples]
    space = Basis(field)
    unit = space.add([int(i == j) for d in sizes for i in range(d) for j in range(d)])
    pending = deque([unit]) if unit is not None else deque()
    while pending:
        vector = pending.popleft()
        for generator in range(len(tuples[0])):
            values = []
            offset = 0
            for w, size in zip(tuples, sizes):
                for a in range(size):
                    for b in range(size):
                        entry = sum(vector[offset+a*size+j]*w[generator][j][b] for j in range(size))
                        values.append(field.reduce(entry))
                offset += size*size
            extra = space.add(values)
            if extra is not None:
                pending.append(extra)
    return space


def included(first, second, field):
    # All simple factors of the second must occur in the first orbit.
    # Its annihilator is the restriction kernel in the joint finite image.
    # A faithful finite module contains all simple types of its image:
    # evaluation at a vector-space basis embeds the regular module in a
    # finite direct sum of copies. Nilpotent action tests containment.
    image = joint_image(first + second, field)
    offset = sum(len(w[0])**2 for w in first)
    kernel = [v for pivot, v in image.rows.items() if pivot >= offset]
    for w in second:
        size = len(w[0])
        entries = [v[offset:offset+size*size] for v in kernel]
        if not ideal_kills_by_iteration(entries, size, field):
            return False
        offset += size*size
    return True


def compute(data):
    if not isinstance(data, dict):
        raise ValueError("the input must be a JSON object")
    p = data["prime"]
    if p is None:
        raise ValueError("supply a prime field")
    field = Field(p)
    op = data["operation"]
    if op == "cubic":
        c = matrix(data["C"], field)
        m, r = len(c), len(c[0])
        d = rank(c, p)
        degree = natural(data.get("degree", 8))
        result = {"rank": d, "degree_one_generators": r+1, "free_degree_two_generators": m-d,
                  "normal_form_relations": [f"x^2 z{i} = z{i} x^2" for i in range(1, d+1)],
                  "residue_field_ext_dimensions": predicted(r, m, d, degree),
                  "profile_poset": "finite-twist localization theorem" if d == r == m else "two elements"}
        if data.get("direct_bar"):
            result["direct_bar_dimensions"] = reduced_ranks(c, p, degree)
        return result
    if op not in ("word", "compare"):
        raise ValueError("operation must be cubic, word or compare")
    n = natural(data["n"], 1)
    q = matrix(data["Q"], field, square=True)
    if rank(q, p) < len(q):
        raise ValueError("Q must be invertible")
    if op == "word":
        word = tuple(natural(a) for a in data["word"])
        if any(a > len(q) for a in word):
            raise ValueError("use 0 for x and i for z_i, starting at 1")
        terms = normal_form_engine(n, q, p)(word)
        return {"normal_form": [{"word": w, "coefficient": c} for w, c in sorted(terms.items())],
                "open_free_rank": n*len(q)}
    orbits = [twist_orbit(read_tuple(data[name], n*len(q), field), n, q, p)
              for name in ("first", "second")]
    ab = included(*orbits, field)
    ba = included(*reversed(orbits), field)
    return {"first_subset_of_second": ab, "second_subset_of_first": ba,
            "equal": ab and ba, "framed_orbit_lengths": [len(w) for w in orbits]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input")
    args = parser.parse_args()
    try:
        if args.input == "-":
            data = json.load(sys.stdin)
        else:
            with open(args.input, encoding="utf-8") as stream:
                data = json.load(stream)
        print(json.dumps(compute(data), indent=2))
    except (ValueError, OSError, KeyError, TypeError, ZeroDivisionError) as error:
        parser.exit(2, f"Invalid input: {error}\n")


if __name__ == "__main__":
    main()
