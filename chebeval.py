import numpy as np


def chebeval(x, p, interval=(-1., 1.), doVelocity=True):
    """!Evaluate a Chebyshev series and first derivative at points x.

    If p is of length n + 1, this function returns:
    y_hat(x) = p_0 * T_0(x*) + p_1 * T_1(x*) + ... + p_n * T_n(x*)
    where T_n(x*) are the orthogonal Chebyshev polynomials of the
    first kind, defined on the interval [-1, 1] and p_n are the
    coefficients. The scaled variable x* is defined on the [-1, 1]
    interval such that (x*) = (2*x - a - b)/(b - a), and x is defined
    on the [a, b] interval.

    \param[in] x: scalar or array of points at which to evaluate the polynomial
    \param[in] p: array-like, chebyshev polynomial coefficients, as returned chebfit.
    \param[in] interval: 2-element list/tuple of bounding the x-interval
                         on which the Chebyshev coefficients were fit.
    \param[in] doVelocity: bool. Do compute the first derivative at points x?
    """
    if len(interval) != 2:
        raise RuntimeError("interval must have length 2")

    intervalBegin = np.float(interval[0])
    intervalEnd = np.float(interval[-1])
    t = 2.*np.array(x, dtype=np.float64) - intervalBegin - intervalEnd
    t /= intervalEnd - intervalBegin

    y = 0.
    v = 0.
    y0 = np.ones_like(t)
    y1 = t
    v0 = np.zeros_like(t)
    v1 = np.ones_like(t)
    v2 = 4.*t
    t = 2.*t
    N = len(p)

    if doVelocity:
        for i in np.arange(0, N, 2):
            if i == N - 1:
                y1 = 0.
                v1 = 0.
            j = min(i + 1, N - 1)

            y += p[i]*y0 + p[j]*y1
            v += p[i]*v0 + p[j]*v1

            y2 = t*y1 - y0
            y3 = t*y2 - y1
            v2 = t*v1 - v0 + 2*y1
            v3 = t*v2 - v1 + 2*y2

            y0 = y2
            y1 = y3
            v0 = v2
            v1 = v3

        return y, 2*v/(intervalEnd - intervalBegin)
    else:
        for i in np.arange(0, N, 2):
            if i == N - 1:
                y1 = 0.
            j = min((i + 1), (N - 1))
            y += p[i]*y0 + p[j]*y1
            y0 = t*y1 - y0
            y1 = t*y0 - y1

        return y, None
