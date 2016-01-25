"""chebgrid.py: Estimate Chebyshev coefficients of a function on a grid

Implementation of Newhall, X. X. 1989, Celestial Mechanics, 45, p. 305-310
"""
import numpy as np
import chebeval as cec
import warnings


def makeChebMatrix(nPoints, nPoly, weight=0.16):
    """Compute C1^(-1)C2 using Newhall89 approach.

    Utility function for fitting chebyshev polynomials to x(t) and dx/dt(t) forcing
    equality at the end points.  This function computes the matrix (C1^(-1)C2).
    Multiplying this matrix by the x and dx/dt values to be fit produces the chebyshev
    coefficient. This function need only be called once for a given polynomial degree and
    number of points.

    The matrices returned are of shape(nPoints+1)x(nPoly).
    The coefficients fitting the nPoints+1 points, X, are found by:
    A = xMultiplier * x  +  dxMultiplier * dxdt if derivative information is known, or
    A = xMultiplier * x  if no derivative information is known.
    The xMultiplier matrices are different, depending on whether derivative information is known.
    Use function makeChebMatrixOnlyX if derviative is not known.
    See Newhall, X. X. 1989, Celestial Mechanics, 45, p. 305-310 for details.

    \param[in] nPoints: Number of point to be fits. Must be greater than 2.
    \param[in] nPoly:   Numper of polynomial terms. Polynomial degree + 1
    \param[in] weight:  Weight to allow control of relative effectos of position and velocity
                        values. Newhall80 found best results are obtained with velocity weighted
                        at 0.4 relative to position, giving W the form (1.0, 0.16, 1.0, 0.16,...)
    \return xMultiplier:  C1^(-1)C2 even rows of shape (nPoints+1)x(nPoly)
                      to be multiplied by x values.
            dxMultiplier: C1^(-1)C2 odd rows of shape (nPoints+1)x(nPoly)
                      to be multiplied by dx/dy values
    """
    tmat = np.zeros([nPoints, nPoly])
    tdot = np.zeros([nPoints, nPoly])

    cj = np.zeros([nPoly])
    xj = np.linspace(1, -1, nPoints)

    for i in np.arange(0, nPoly):
        cj[:] = 0
        cj[i] = 1
        y, v = cec.chebeval(xj, cj)
        tmat[:, i] = y
        tdot[:, i] = v

    # make matrix T*W
    tw = np.zeros([nPoly, nPoints, 2])
    tw[:, :, 0] = tmat.transpose()
    tw[:, :, 1] = tdot.transpose()*weight

    # make matrix T*WT
    twt = np.dot(tw[:, :, 0], tmat) + np.dot(tw[:, :, 1], tdot)
    tw = tw.reshape(nPoly, 2*nPoints)

    # insert matrix T*W in matrix C2
    c2 = np.zeros([nPoly + 4, 2 * nPoints])
    c2[0:nPoly] = tw
    c2[nPoly, 0] = 1
    c2[nPoly + 1, 1] = 1
    c2[nPoly + 2, -2] = 1
    c2[nPoly + 3, -1] = 1

    # insert matrix T*WT in matrix C1
    c1 = np.zeros([nPoly + 4, nPoly + 4])
    c1[0:nPoly, 0:nPoly] = twt
    c1[nPoly + 0, 0:nPoly] = tmat[0]
    c1[nPoly + 1, 0:nPoly] = tdot[0]
    c1[nPoly + 2, 0:nPoly] = tmat[-1]
    c1[nPoly + 3, 0:nPoly] = tdot[-1]

    c1[0:nPoly, nPoly:] = c1[nPoly:, 0:nPoly].transpose()

    c1inv = np.linalg.inv(c1)
    c1c2 = np.dot(c1inv, c2)
    c1c2 = c1c2.reshape(nPoly + 4, nPoints, 2)
    c1c2 = c1c2[:, ::-1, :]
    c1c2 = c1c2.reshape(nPoly + 4, 2*nPoints)

    # separate even rows for x, and odd rows for dx/dt
    return c1c2[0:nPoly, 0::2], c1c2[0:nPoly, 1::2]


def makeChebMatrixOnlyX(nPoints, nPoly):
    """Compute C1^(-1)C2 using Newhall89 approach without dx/dt

    Compute xMultiplier using only the equality constraint of the x-values at the endpoints.
    To be used when first derivatives are not available.
    If chebyshev approximations are strung together piecewise only the x-values
    and not the first derivatives will be continuous at the boundaries.
    Multiplying this matrix by the x-values to be fit produces the chebyshev
    coefficients. This function need only be called once for a given polynomial degree and
    number of points. See Newhall, X. X. 1989, Celestial Mechanics, 45, p. 305-310.

    @param [in] nPoints: Number of point to be fits. Must be greater than 2.
    @param [in] nPoly:   Numper of polynomial terms. Polynomial degree + 1
    \return xMultiplier: Even rows of C1^(-1)C2 w/ shape (nPoints+1)x(nPoly)
                         to be multiplied by x values
    """

    tmat = np.zeros([nPoints, nPoly])
    cj = np.zeros([nPoly])
    xj = np.linspace(1, -1, nPoints)
    for i in range(0, nPoly):
        cj[:] = 0
        cj[i] = 1
        tmat[:, i], v = cec.chebeval(xj, cj)

    # Augment matrix T to get matrix C2
    c2 = np.zeros([nPoly + 2, nPoints])
    c2[0:nPoly] = tmat.transpose()
    c2[nPoly, 0] = 1
    c2[nPoly + 1, nPoints - 1] = 1

    # Augment matrix T*WT to get the matrix C1
    c1 = np.zeros([nPoly + 2, nPoly + 2])
    c1[0:nPoly, 0:nPoly] = np.dot(tmat.transpose(), tmat)
    c1[nPoly + 0, 0:nPoly] = tmat[0]
    c1[nPoly + 1, 0:nPoly] = tmat[-1]
    c1[0:nPoly, nPoly:] = c1[nPoly:, 0:nPoly].transpose()

    c1inv = np.linalg.inv(c1)
    # C1^(-1) C2
    c1c2 = np.dot(c1inv, c2)

    c1c2 = c1c2.reshape(nPoly + 2, nPoints)
    c1c2 = c1c2[:, ::-1]
    return c1c2[0:nPoly]


def chebfit(t, x, dxdt=None, xMultiplier=None, dxMultiplier=None, nPoly=7):
    """!Fit Chebyshev polynomial constrained at endpoints using Newhall89 approach.

    Return Chebyshev coefficients and statistics from fit
    to array of positions (x) and optional velocities (dx/dt).
    If both the function and its derivative are specified, then the value and
    derivative of the interpolating polynomial at the
    endpoints will be exactly equal to the input endpoint values.
    Many approximations may be piecewise strung together and the function value
    and its first derivative will be continuous across boundaries. If derivatives
    are not provided, only the function value will be continuous across boundaries.

    If xMultiplier and dxMultiplier are not provided or
    are an inappropriate shape for t and x, they will be recomputed.
    See Newhall, X. X. 1989, Celestial Mechanics, 45, p. 305-310
    for details.

    \param[in] t:      Array of regularly sampled independent variable (e.g. time)
    \param[in] x:      Array of regularly sampled dependent variable (e.g. declination)
    \param[in] dxdt:   Optionally, array of first derivatives of x with respect to t,
                       at the same grid points. (e.g. sky velocity ddecl/dt)
    \param[in] xMultiplier:  Optional 2D Matrix with rows of C1^(-1)C2 corresponding to x.
                       Use makeChebMatrix to compute
    \param[in] dxMultiplier: Optional 2D Matrix with rows of C1^(-1)C2 corresponding to dx/dt.
                       Use makeChebMatrix to compute
    \param[in] nPoly:  Number of polynomial terms. Degree + 1.  Must be >=2 and < 2*(npoints + 1),
                       when derivative information is specified or less than npoints + 1, when no
                       derivative information is specified.

    \return: tuple containing
                Array of chebyshev coefficients with length=nPoly,
                Array of residuals of the tabulated function x minus the approximated function
                float: rms residuals
    """
    nPoints = len(t)
    if len(x) != nPoints:
        raise ValueError("length of x (%s) != length of t (%s)" % (len(x), nPoints))
    if dxdt is None:
        if nPoly >= nPoints:
            raise RuntimeError('Without velocity constraints, nPoly must be less than %s' % (nPoints + 1))
        if nPoly < 2:
            raise RuntimeError('Without velocity constraints, nPoly must be greater than 2')
    else:
        if nPoly >= 2 * (nPoints + 1):
            raise RuntimeError('nPoly must be less than %s' % (2*(nPoints + 1)))
        if nPoly < 4:
            raise RuntimeError('nPoly must be greater than 4')

    # Recompute C1invX2 if xMultiplier and dxMultiplier are None or
    # they are not appropriate for sizes of input positions and velocities.

    if xMultiplier is None:
        redoX = True
    else:
        redoX = (xMultiplier.shape[1] != nPoints) | (xMultiplier.shape[0] != nPoly)

    if dxMultiplier is None:
        redoV = True
    else:
        redoV = (dxMultiplier.shape[1] != nPoints) | (dxMultiplier.shape[0] != nPoly)

    if (dxdt is None) & redoX:
        warnings.warn("Chebfit making x-only multiplier for npoints=%s, coeff=%s." % (nPoints, nPoly) +
                      "Routine could be sped up by precumputing matrix using makeChebMatrixOnlyX()")
        xMultiplier = makeChebMatrixOnlyX(nPoints, nPoly)

    if (dxdt is not None) & (redoV | redoX):
        warnings.warn("Chebfit making multipliers for npoints=%s, coeff=%s." % (nPoints, nPoly) +
                      "Routine could be sped up by precumputing matrix using makeChebMatrix()")
        xMultiplier, dxMultiplier = makeChebMatrix(nPoints, nPoly)

    if x.size != nPoints:
        raise RuntimeError("Not enough elements in X")

    tInterval = np.array([t[0], t[-1]]) - t[0]
    tScaled = t - t[0]

    # Compute the X portion of the coefficients
    a_n = np.dot(xMultiplier, x)

    # Compute statistics
    # for x and dxdt if it is available
    if dxdt is not None:
        a_n = a_n + np.dot(dxMultiplier, dxdt*(tInterval[1] - tInterval[0])/2.)
        xApprox, dxApprox = cec.chebeval(tScaled, a_n, interval=tInterval)
        dresiduals = dxdt - dxApprox
        dse = np.sum(dresiduals**2)
        drms = np.sqrt(dse/(nPoints - 1))
    else:
        # Statistics for x only
        xApprox, _ = cec.chebeval(tScaled, a_n, interval=tInterval, doVelocity=False)

    residuals = x - xApprox
    se = np.sum(residuals**2)
    rms = np.sqrt(se/(nPoints - 1))

    return a_n, residuals, rms
