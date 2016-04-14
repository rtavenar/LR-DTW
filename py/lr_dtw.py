from scipy.spatial.distance import cdist
import numpy

__author__ = 'Romain Tavenard romain.tavenard[at]univ-rennes2.fr'

UP = 0
RIGHT = 1
DIAGONAL = 2


# @deprecated (+buggy?) function
def get_probas(cost_up, cost_right, cost_diagonal, gamma):
    probas = numpy.zeros((3,)) * numpy.nan
    vec_p1 = numpy.arange(0., 1.01, .01)
    vec_p2 = numpy.arange(0., 1.01, .01)
    p1, p2 = numpy.meshgrid(vec_p1, vec_p2)
    C = cost_up * p1 + cost_right * p2 + cost_diagonal * (1 - p1 - p2) + gamma * (
    p1 ** 2 + p2 ** 2 + (1 - p1 - p2) ** 2)
    C[p1 < 0] = C[p2 < 0] = C[p1 + p2 > 1] = numpy.inf
    Cmin = numpy.min(C)
    p1_best = p1[C == Cmin][0]
    p2_best = p2[C == Cmin][0]
    probas[UP] = p1_best
    probas[DIAGONAL] = p2_best
    probas[RIGHT] = 1 - p1_best - p2_best
    assert numpy.alltrue(numpy.logical_and(probas >= -1e-10, probas <= 1 + 1e-10)), probas
    assert numpy.sum(probas) <= 1 + 1e-10, probas
    return probas


def get_probas_formula(cost_up, cost_right, cost_diagonal, gamma, entropy_penalized=False):
    if gamma < 1e-12:
        probas = numpy.zeros((3,))
        min_val = min(cost_up, cost_right, cost_diagonal)
        if cost_up == min_val:
            probas[UP] = 1.
        if cost_diagonal == min_val:
            probas[DIAGONAL] = 1.
        if cost_right == min_val:
            probas[RIGHT] = 1.
        return probas / numpy.sum(probas)
    probas = numpy.zeros((3,)) * numpy.nan
    if entropy_penalized:
        p_up = 1. / (1. + numpy.power(2., (cost_up - cost_right) / gamma) +
                     numpy.power(2., (cost_up - cost_diagonal) / gamma))
        p_right = 1. / (1. + numpy.power(2., (cost_right - cost_up) / gamma) +
                        numpy.power(2., (cost_right - cost_diagonal) / gamma))
    else:
        p_up = 1. / 3. * (1. + (cost_right + cost_diagonal - 2 * cost_up) / (2 * gamma))
        p_right = 1. / 3. * (1. + (cost_up + cost_diagonal - 2 * cost_right) / (2 * gamma))
    if p_up >= 0. and p_right >= 0. and p_up + p_right <= 1.:
        probas[UP], probas[RIGHT], probas[DIAGONAL] = p_up, p_right, 1 - p_up - p_right
    elif p_up < 0:
        p_up = 0.
        if entropy_penalized:
            p_right = 1. / (1. + numpy.power(2., (cost_right - cost_diagonal) / gamma))
        else:
            p_right = .5 * (1. + (cost_diagonal - cost_right) / (2 * gamma))
        if 0 <= p_right <= 1:
            probas[UP], probas[RIGHT], probas[DIAGONAL] = p_up, p_right, 1 - p_up - p_right
        elif p_right < 0.:
            probas[UP], probas[RIGHT], probas[DIAGONAL] = 0., 0., 1.
        else:
            probas[UP], probas[RIGHT], probas[DIAGONAL] = 0., 1., 0.
    elif p_right < 0:
        p_right = 0.
        if entropy_penalized:
            p_up = 1. / (1. + numpy.power(2., (cost_up - cost_diagonal) / gamma))
        else:
            p_up = .5 * (1. + (cost_diagonal - cost_up) / (2 * gamma))
        if 0 <= p_up <= 1:
            probas[UP], probas[RIGHT], probas[DIAGONAL] = p_up, p_right, 1 - p_up - p_right
        elif p_up < 0.:
            probas[UP], probas[RIGHT], probas[DIAGONAL] = 0., 0., 1.
        else:
            probas[UP], probas[RIGHT], probas[DIAGONAL] = 1., 0., 0.
    else:
        p_diagonal = 0.
        if entropy_penalized:
            p_up = 1. / (1. + numpy.power(2., (cost_up - cost_right) / gamma))
        else:
            p_up = .5 * (1. + (cost_right - cost_up) / (2 * gamma))
        if 0 <= p_up <= 1:
            probas[UP], probas[RIGHT], probas[DIAGONAL] = p_up, 1 - p_diagonal - p_up, p_diagonal
        elif p_up < 0.:
            probas[UP], probas[RIGHT], probas[DIAGONAL] = 0., 1., 0.
        else:
            probas[UP], probas[RIGHT], probas[DIAGONAL] = 1., 0., 0.
    assert numpy.alltrue(numpy.logical_and(probas >= -1e-10, probas <= 1 + 1e-10)), probas
    assert numpy.sum(probas) <= 1 + 1e-10, probas
    return probas


def lr_dtw(s_x=None, s_y=None, mat_dist=None, gamma=0., entropy_penalized=False):
    if mat_dist is None:
        assert s_x is not None and s_y is not None, "If mat_dist is not given, both time series should be provided"
        mat_dist = cdist(s_x, s_y)
    n, m = mat_dist.shape
    mat_cost = numpy.zeros((n, m))
    probas = numpy.zeros((n, m, 3))
    for i in range(1, n):
        probas[i, 0][UP] = 1.
        mat_cost[i, 0] = mat_cost[i - 1, 0] + mat_dist[i, 0]
    for j in range(1, m):
        probas[0, j][RIGHT] = 1.
        mat_cost[0, j] = mat_cost[0, j - 1] + mat_dist[0, j]
    for i in range(1, n):
        for j in range(1, m):
            probas[i, j] = get_probas_formula(cost_up=mat_cost[i - 1, j], cost_right=mat_cost[i, j - 1],
                                              cost_diagonal=mat_cost[i - 1, j - 1], gamma=gamma,
                                              entropy_penalized=entropy_penalized)
            mat_cost[i, j] = probas[i, j][UP] * mat_cost[i - 1, j] + probas[i, j][RIGHT] * mat_cost[i, j - 1] + \
                             probas[i, j][DIAGONAL] * mat_cost[i - 1, j - 1] + mat_dist[i, j]
    return mat_cost[-1, -1], probas


def lr_dtw_backtrace(probas):
    n, m = probas.shape[:2]
    mat_probas = numpy.zeros((n, m))
    mat_probas[-1, -1] = 1.
    for i in range(n - 2, -1, -1):
        mat_probas[i, -1] = mat_probas[i + 1, -1] * probas[i + 1, -1][UP]
    for j in range(m - 2, -1, -1):
        mat_probas[-1, j] = mat_probas[-1, j + 1] * probas[-1, j + 1][RIGHT]
    for i in range(n - 2, -1, -1):
        for j in range(m - 2, -1, -1):
            mat_probas[i, j] = mat_probas[i + 1, j] * probas[i + 1, j][UP] + \
                               mat_probas[i, j + 1] * probas[i, j + 1][RIGHT] + \
                               mat_probas[i + 1, j + 1] * probas[i + 1, j + 1][DIAGONAL]
    return mat_probas


if __name__ == "__main__":
    for gamma in [0., 1., 10., 100., 1000.]:
        print(gamma, get_probas(10, 15, 20, gamma))
