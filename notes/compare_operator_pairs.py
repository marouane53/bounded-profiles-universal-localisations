#!/usr/bin/env python3
"""Exact inclusion and equality for the operator-pair bounded profiles.

Pass a JSON file (or - for stdin) with keys first and second, each [T0,T1],
and optional prime. Without prime the field is Q. Fractions can be strings.
The two pairs may have different dimensions.
"""
import argparse
import json
from collections import deque
from fractions import Fraction
from math import isqrt
import sys


class Field:
    def __init__(self, prime=None):
        if prime is not None:
            if (type(prime) is not int or prime < 2
                    or any(prime % d == 0 for d in range(2, isqrt(prime) + 1))):
                raise ValueError('prime must be a prime integer')
        self.prime = prime

    def parse(self, value):
        if type(value) not in (str, int):
            raise ValueError('use integers or quoted rational numbers')
        value = Fraction(value)
        if self.prime is None:
            return value
        return (value.numerator * pow(value.denominator, -1, self.prime)) % self.prime

    def reduce(self, value):
        return value if self.prime is None else value % self.prime

    def inverse(self, value):
        return Fraction(1, value) if self.prime is None else pow(value, -1, self.prime)


class Basis:
    def __init__(self, field):
        self.field = field
        self.rows = {}

    def add(self, vector):
        """Return the new normalized basis vector, or None when dependent."""
        row = list(vector)
        for i in range(len(row)):
            a = row[i]
            if not a:
                continue
            if i not in self.rows:
                factor = self.field.inverse(a)
                row = tuple(self.field.reduce(v * factor) for v in row)
                self.rows[i] = row
                return row
            old = self.rows[i]
            for j in range(i, len(row)):
                row[j] = self.field.reduce(row[j] - a * old[j])
        return None


def read_pair(pair, field):
    if not isinstance(pair, list) or len(pair) != 2 or not isinstance(pair[0], list):
        raise ValueError('a pair must be a list of two square matrices')
    size = len(pair[0])
    for matrix in pair:
        if (not isinstance(matrix, list) or len(matrix) != size
                or any(not isinstance(row, list) or len(row) != size for row in matrix)):
            raise ValueError('the two matrices in a pair must be square and of the same size')
    return tuple(tuple(tuple(field.parse(x) for x in row) for row in m) for m in pair)


def joint_image(summands, field):
    sizes = [len(pair[0]) for pair in summands]
    identity = [int(i == j) for d in sizes for i in range(d) for j in range(d)]
    space = Basis(field)
    unit = space.add(identity)
    todo = deque([unit]) if unit is not None else deque()
    while todo:
        row = todo.popleft()
        for generator in range(2):
            offset = 0
            prod = []
            for pair, d in zip(summands, sizes):
                matrix = pair[generator]
                for i in range(d):
                    for j in range(d):
                        val = sum(row[offset+i*d+k]*matrix[k][j] for k in range(d))
                        prod.append(field.reduce(val))
                offset += d*d
            extra = space.add(prod)
            if extra is not None:
                todo.append(extra)
    return space


def ideal_kills_by_iteration(ideal, size, field):
    current = [tuple(int(i == j) for j in range(size)) for i in range(size)]
    while current:
        new = Basis(field)
        for m in ideal:
            for v in current:
                new.add([field.reduce(sum(m[i*size+j]*v[j] for j in range(size)))
                         for i in range(size)])
        if len(new.rows) == len(current):
            return False
        current = list(new.rows.values())
    return True


def profile_included(first, second, field):
    summands = [first, first[::-1], second, second[::-1]]
    image = joint_image(summands, field)
    cutoff = 2*len(first[0])**2
    size = len(second[0])
    # Leading pivots after the first two summands span the restriction kernel.
    kernel = [v for pivot, v in image.rows.items() if pivot >= cutoff]
    for offset in [cutoff, cutoff + size*size]:
        ideal = [v[offset:offset+size*size] for v in kernel]
        if not ideal_kills_by_iteration(ideal, size, field):
            return False
    return True


def compare(data):
    field = Field(data.get('prime'))
    a, b = (read_pair(data[name], field) for name in ('first', 'second'))
    ab = profile_included(a, b, field)
    ba = profile_included(b, a, field)
    return {'first_subset_of_second': ab, 'second_subset_of_first': ba, 'equal': ab and ba}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', help='JSON file or - for stdin')
    args = parser.parse_args()
    try:
        if args.input == '-':
            data = json.load(sys.stdin)
        else:
            with open(args.input, encoding='utf-8') as stream:
                data = json.load(stream)
        if not isinstance(data, dict):
            raise ValueError('the JSON input must be an object')
        print(json.dumps(compare(data), indent=2))
    except (OSError, KeyError, ValueError, ZeroDivisionError) as error:
        parser.exit(2, f'Invalid input: {error}\n')


if __name__ == '__main__':
    main()
