"""The surface, in one place.

Every certificate in this repository derives the defining data of the surface
from this module rather than restating it, so that a change here would turn
the checks red instead of silently disagreeing with them.

    X = { Z in P^5 : F_0(Z) = F_1(Z) = F_2(Z) = 0 },
    F_k(Z) = sum_{j=0}^{5} mu_j^k z_j^2,     k = 0, 1, 2.

The six parameters are pairwise distinct positive integers. Distinctness is
what makes every 3x3 Vandermonde block V_S invertible and every 2x2 block
nonzero, and those two facts are exactly what the smoothness and coverage
arguments consume.

The equations are LINEAR in the squares w_j = z_j^2. Writing V_full for the
3x6 matrix (mu_j^k), the surface imposes V_full . w = 0; resolving a triple S
of coordinates against the remaining three is the chart construction used
throughout.
"""

MU_INT = (1, 2, 3, 5, 7, 11)


def quadric(k: int, w):
    """F_k evaluated on the squares w = (w_0, ..., w_5), in exact arithmetic."""
    return sum((mu ** k) * wj for mu, wj in zip(MU_INT, w))


def vandermonde(triple):
    """The 3x3 block V_S = (mu_s^k) for k = 0,1,2 and s in the triple."""
    return [[MU_INT[s] ** k for s in triple] for k in range(3)]
