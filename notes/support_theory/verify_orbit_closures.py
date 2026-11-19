"""Orbit/return computations independent of the geometric proof. Standard library."""
from itertools import product

from check_lib import require_asserts


def polynomial_difference(a, b, prime):
    return [(x - y) % prime for x, y in zip(a, b)]


def frobenius_orbit(prime, steps):
    # The two slots for each generator start equal. Substitution sends
    # (u0,u1,v0,v1) to (u1,t*u0,v1,(t+1)*v0) in F_p[t].
    one = [1] + [0] * (steps + 1)
    u0 = u1 = v0 = v1 = one
    returns = []
    for i in range(2 * steps + 2):
        value = polynomial_difference(polynomial_difference(v0, u0, prime), one, prime)
        if not any(value):
            returns.append(i)
        times_t = [0] + u0[:-1]
        times_t_plus_one = [(a + b) % prime for a, b in zip(v0, [0] + v0[:-1])]
        u0, u1, v0, v1 = u1, times_t, v1, times_t_plus_one

    powers = []
    j = 1
    while j <= steps:
        powers.extend((2 * j, 2 * j + 1))
        j *= prime
    assert returns == powers
    return len(returns)


def det(a):
    return a[0] * a[3] - a[1] * a[2]


def multiply(a, b):
    return tuple(sum(a[2*i+k] * b[2*k+j] for k in (0, 1))
                 for i in (0, 1) for j in (0, 1))


def translate(a, b, c):
    return tuple(x + c*y for x, y in zip(a, b))


def nonlinear_return_test():
    # V -> V+U^2 is a nonlinear free-algebra automorphism.
    # Work with integral 2x2 framed tuples in both time directions.
    matrices = list(product(range(2), repeat=4))
    total = 0
    for u, v in product(matrices, repeat=2):
        step = multiply(u, u)
        start = translate(v, u, 1)
        fp, fm, f0 = det(translate(start, step, 1)), det(translate(start, step, -1)), det(start)
        linear = (fp - fm) // 2
        quadratic = (fp + fm - 2*f0) // 2
        current = translate(v, step, -12)
        failures = []
        for time in range(-12, 13):
            value = det(translate(current, u, 1))
            assert value == f0 + time*linear + time*time*quadratic
            if value == 0:
                failures.append(time)
            current = translate(current, step, 1)
        identically_zero = (f0, linear, quadratic) == (0, 0, 0)
        assert len(failures) == 25 if identically_zero else len(failures) <= 2
        total += 1
    # Two distinct components, only one of which is contained in the bad divisor.
    u, v = 1, 0
    hits = []
    for i in range(50):
        if u == 1:
            hits.append(v)
        u, v = -u, v + 1
    assert hits == list(range(0, 50, 2))
    return total


def fixed_dimension_separators():
    # Simple tuples (c*S,T) with diagonal T and nonzero off-diagonals in S.
    # Scaling by 2 gives dense orbits in each line; 3 is not a power of 2.
    # V-2 and V-3 distinguish the two conjugacy-saturated varieties.
    s = (0, 1, 1, 0)
    identity = (1, 0, 0, 1)
    for eigenvalue in (2, 3):
        v = (1, 0, 0, eigenvalue)
        for c in (1, 3):
            for j in range(20):
                u = tuple(c*(2**j)*a for a in s)
                commutator = translate(multiply(u, v), multiply(v, u), -1)
                assert det(commutator) != 0
                for root in (2, 3):
                    assert (det(translate(v, identity, -root)) == 0) == (root == eigenvalue)


if __name__ == "__main__":
    require_asserts()
    counts = {p: frobenius_orbit(p, 256) for p in (2, 3, 5, 7, 11)}
    print("Polynomial orbits through degree 513 over F_p(t). Counts of failed twists:", counts)
    print(f"{nonlinear_return_test()} integral matrix pairs under a nonlinear automorphism, "
          "25 twists each; also checked a reducible orbit with alternating returns.")
    fixed_dimension_separators()
    print("Checked conjugacy-invariant polynomial witnesses for two-dimensional simple tuples.")
