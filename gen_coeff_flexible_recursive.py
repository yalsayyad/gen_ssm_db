"""
Compute ephemerides with openorb and output their
Chebyshev polynomial coefficients

Usage:
python gen_coeff_flexible_recursive.py <filename> <start_time_tai> <days> <coefficients> <totaldays>

Examples: python gen_coeff_flexible_recursive.py tests/S1_00.zzz_49353.des 49343 30 14 30

remember to setup the library environment variables first:
>export OORB_CONF=/share/pogo3/krughoff/shared/oorb/main/oorb.conf
>export OORB_DATA=/share/pogo3/krughoff/shared/oorb/data/
"""

import numpy as np
import chebfit as cg
import chebeval as ce
import sys
import movingObject as mo
import math
import pyoorb as oo
DEBUG = False


def sphericalDistance(origin, destination):
    lat1, lon1 = origin
    lat2, lon2 = destination
    radius = 1.0
    dlat = math.radians(lat2-lat1)
    dlon = math.radians(lon2-lon1)
    a = math.sin(dlat/2) * math.sin(dlat/2) + math.cos(math.radians(lat1)) \
        * math.cos(math.radians(lat2)) * math.sin(dlon/2) * math.sin(dlon/2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    d = radius * c
    return d*180.0/np.pi


def get_coeffs_position(t, ra, dec, dracoorddt, ddecdt, ng, npo, coeff, multiplier=(None, None)):
    p, dec_resid, dec_rms = cg.chebfit(t, dec, ddecdt,
                                       xMultiplier=multiplier[0], dxMultiplier=multiplier[1], nPoly=coeff)
    rap, ra_resid, ra_rms = cg.chebfit(t, ra, dracoorddt,
                                       xMultiplier=multiplier[0], dxMultiplier=multiplier[1], nPoly=coeff)
    ra_real_resid = (ra_resid)*np.cos(np.pi*dec/180.)
    return p, rap, 3600.0*1000.0*np.max(np.sqrt(dec_resid**2 + ra_real_resid**2))


def get_coeffs_vmag(t, vmag, ng, npo, coeff, multiplier=None):
    p, resid, rms = cg.chebfit(t, vmag, None, xMultiplier=multiplier, nPoly=coeff)
    return p, np.max(np.abs(resid))


def get_coeffs_se(t, se, ng, npo, coeff, multiplier=None):
    p, resid, rms = cg.chebfit(t, se, None, xMultiplier=multiplier, nPoly=coeff)
    return p, np.max(np.abs(resid))


def get_coeffs_dist(t, dist, distdt, ng, npo, coeff, multiplier=(None, None)):
    p, resid, rms = cg.chebfit(t, dist, distdt,
                               xMultiplier=multiplier[0], dxMultiplier=multiplier[1], nPoly=coeff)
    return p, np.max(np.abs(resid))


def get_coeffs_dist2(t, dist, ng, npo, coeff, multiplier=None):
    p, resid, rms = cg.chebfit(t, dist, None, xMultiplier=multiplier, nPoly=coeff)
    return p, np.max(np.abs(resid))


def get_granularity(distance):
    """
    if distance is 0.8 degrees/day  treat same as MBA
                 < 1.6 degrees/day  try gen 1 day at 64 points per day.
                 < 3.2 deg/day      try gen 0.5 day at 128 points per day
                 <6.4 deg/day       try gen 0.25 day at 256 points per day
                 < 12.8 deg.day     try gen 0.125 day at 512 points per day
                 < 25.6 deg/day     try gen 0.0625 day at 1024 points per day
                 < 51.2 deg/day     try gen 0.03125 day at 2048 points per day
                                    try gen 0.015625 day at 4096 points per day
    ngran = 64 #always ngran = int(range/timestep)
    """
    if distance < 0.8:
        timestep = 0.03125  # 1/32 day
        length = 2.         # days 2
    elif distance < 1.6:
        timestep = 0.015625  # 1/64 day
        length = 1.          # days 1
    elif distance < 3.2:
        timestep = 0.0078125  # 1/128 day
        length = 0.5          # days 1/2
    elif distance < 6.4:
        timestep = 0.00390625  # 1/256 day
        length = 0.25          # days 1/4
    elif distance < 12.8:
        timestep = 0.001953125  # 1/512 day
        length = 0.125          # days 1/8
    elif distance < 25.6:
        timestep = 0.0009765625  # 1/1024 day
        length = 0.0625          # days 1/16
    elif distance < 51.2:
        timestep = 0.00048828125  # 1/2048 day
        length = 0.03125          # days 1/32
    elif distance < 102.4:
        timestep = 0.000244140625  # 1/4096 day
        length = 0.015625          # days 1/64
    else:  # fastest it can go
        timestep = 0.0001220703125  # 1/8192 day
        length = 0.0078125          # days  1/128
    return timestep, length, 64


def three_sixy_to_neg(element, min, max):
    if (min < 100) & (max > 270):
        if element > 270:
            return element - 360.
        else:
            return element
    else:
        return element

v_360_to_neg = np.vectorize(three_sixy_to_neg)


def getEphem(movobj, startMjd, length, dt):
    mjdTaiList = list(np.arange(startMjd, startMjd + length + dt, dt))
    t = np.empty(len(mjdTaiList))
    ra = np.empty(len(mjdTaiList))
    dec = np.empty(len(mjdTaiList))
    dracoorddt = np.empty(len(mjdTaiList))
    ddecdt = np.empty(len(mjdTaiList))
    vmag = np.empty(len(mjdTaiList))
    dist = np.empty(len(mjdTaiList))
    se = np.empty(len(mjdTaiList))
    movobj.calcEphemeris(mjdTaiList, obscode=807, eph_timescale=4.0)
    # if DEBUG==True: print  "end calc ephemeris", time.time() - starttime
    for j, mjdTai in enumerate(mjdTaiList):
        t[j], ra[j], dec[j], dracoorddt[j], ddecdt[j], vmag[j], dist[j], se[j] = \
            movobj.Ephemerides[movobj.mjdTaiStr(mjdTai)].getPosition()
    return t, ra, dec, dracoorddt, ddecdt, vmag, dist, se


def plotRaDecT(ra, dec, t, cra, cdec):
    import matplotlib.pyplot as plt
    trange = np.arange(t[0], t[-1], 0.01)
    raEval, vraEval = ce.chebeval(trange, cra, interval=np.array([t[0], t[-1]]))
    plt.plot(t, ra, 'k.')
    plt.plot(trange, raEval, 'r-')
    plt.ylabel('RA')
    plt.show()
    decEval, vdecEval = ce.chebeval(trange, cdec, interval=np.array([t[0], t[-1]]))
    plt.plot(t, dec, 'k.')
    plt.plot(trange, decEval, 'b-')
    plt.ylabel('Dec')
    plt.show()


def adjustTimeLengthExtreme(p_resid, dec, timestep, length):
    if p_resid > 1000:
        timestep = timestep/16
        length = length/16
    elif p_resid > 100:
        timestep = timestep/8
        length = length/8
    elif p_resid > 15:
        timestep = timestep/4
        length = length/4
    elif p_resid > 5:
        timestep = timestep/2
        length = length/2
    elif p_resid > 2:
        timestep = timestep/1
        length = length/1
    elif p_resid > 1:
        timestep = timestep/1
        length = length/1
    # cut it in half once more if chance to go over poles
    if dec < -75. or dec > 75.:
        timestep = timestep/2
        length = length/2
    return timestep, length


def adjustTimeLengthNormal(p_resid, dec, timestep, length):
    if p_resid > 5:
        timestep = timestep/4
        length = length/4
    elif p_resid > 2:
        timestep = timestep/2
        length = length/2
    # cut it in half once more if chance to go over poles
    if dec < -75. or dec > 75.:
        timestep = timestep/2
        length = length/2
    return timestep, length


def getGran(ssmid, mymo, start_time, days, coeff):
    # first, lets generate a day and see how far it goes:
    halfofDaysToRun = start_time + int(days/2)
    t, ra, decl, dracoorddt, ddecdt, magV, dist, se = getEphem(mymo, halfofDaysToRun, 1.0, 1.0)
    distance = sphericalDistance([ra[0], decl[0]], [ra[1], decl[1]])
    timestep, length, ngran = get_granularity(distance)
    t, ra, dec, dracoorddt, ddecdt, vmag, dist, se = getEphem(mymo, halfofDaysToRun, length, timestep)
    pdec, pra, p_resid = get_coeffs_position(t, v_360_to_neg(ra, np.min(ra), np.max(ra)), dec,
                                             dracoorddt, ddecdt, 64, 64, coeff)
    if DEBUG:
        print "test presid MAIN is", p_resid
    timestep, length = adjustTimeLengthNormal(p_resid, dec[0], timestep, length)
    return timestep, length, ngran


def breakItDown(ssmid, mymo, start_time, timestep,  days,  ngran,  coeff,  multiplier, CoeffFile,
                ResidualSumfile, Failedfile,  inputfilename):
    timestep = timestep/2.
    length = days/2
    halfofDaysToRun = start_time
    t, ra, dec, dracoorddt, ddecdt, vmag, dist, se = getEphem(mymo, halfofDaysToRun, length, timestep)

    pdec, pra, p_resid = get_coeffs_position(t, v_360_to_neg(ra, np.min(ra), np.max(ra)),
                                             dec, dracoorddt, ddecdt, 64, 64, coeff)
    if DEBUG:
        print "test presid RECURSIVE is", p_resid
    timestep, length = adjustTimeLengthExtreme(p_resid, dec[0], timestep, length)
    doOneRecursiveSegment(ssmid, mymo, start_time, days, coeff, timestep, length, ngran,
                          multiplier, CoeffFile, ResidualSumfile, Failedfile, inputfilename)


def doOneRecursiveSegment(ssmid, mymo, start_time, days, coeff, timestep, length, ngran,
                          multiplier, CoeffFile, ResidualSumfile, Failedfile, inputfilename):
    t, ra, dec, dracoorddt, ddecdt, vmag, dist, se = getEphem(mymo, start_time, days, timestep)
    # now we have our arrays and fit exactly like before
    day0 = 0
    day1 = 64
    rows = 1
    npoint = 64
    while day0 < len(t) - 64:
        pdec, pra, p_resid = get_coeffs_position(t[day0:day1+1],
                                                 v_360_to_neg(ra[day0:day1+1],
                                                              np.min(ra[day0:day1+1]),
                                                              np.max(ra[day0:day1+1])),
                                                 dec[day0:day1+1], dracoorddt[day0:day1+1], ddecdt[day0:day1+1],
                                                 ngran, npoint, coeff, multiplier['POSITION'])
        if p_resid > 2.5:
            if DEBUG:
                print "oh no!:", p_resid
            breakItDown(ssmid, mymo, t[day0], timestep, length, ngran, coeff, multiplier,
                        CoeffFile, ResidualSumfile, Failedfile,  inputfilename)
        else:  # we're good. Print it out to file
            d, d_resid = get_coeffs_dist2(t[day0:day1+1], dist[day0:day1+1], ngran, npoint, 5, multiplier['DIST_X'])
            v, v_resid = get_coeffs_vmag(t[day0:day1+1], vmag[day0:day1+1], ngran, npoint, 9, multiplier['VMAG_X'])
            s, s_resid = get_coeffs_se(t[day0:day1+1], se[day0:day1+1], ngran, npoint, 6, multiplier['SE_X'])
            if np.isnan(v_resid) | np.isnan(d_resid) | np.isnan(s_resid):
                print 'do not print!!'
                print >>Failedfile, "%s %i %.14f %.14f %.14f %.14e %.14e %.14e %.14e %s" % (
                                    ssmid, rows, t[day0], t[day1], t[day1] - t[day0], p_resid, d_resid,
                                    v_resid, s_resid, inputfilename)
            else:
                print >>ResidualSumfile, "%s %i %.14f %.14f %.14f %.14e %.14e %.14e %.14e %s" % (
                    ssmid, rows, t[day0], t[day1], t[day1] - t[day0], p_resid, d_resid,
                    v_resid, s_resid, inputfilename)

                print >>CoeffFile, "%i %s %.10f %.10f %s %s %s %s %s"%(0,
                                    ssmid, t[day0], t[day1],  " ".join('%.14e'%j for j in pra), " ".join('%.14e'%j for j in pdec), 
                                    " ".join('%.7e'%j for j in d), " ".join('%.7e'%j for j in v),   " ".join('%.7e'%j for j in s))

        # advance to the next day if less than 6 points left in month
        day0 = day1
        day1 = day1 + 64
        rows = rows + 1


def doOneMonth(ssmid, mymo, start_time, days, coeff, multiplier, CoeffFile, ResidualSumfile,
               Failedfile, inputfilename, knownGran=True):
    if not knownGran:
        timestep, length, ngran = getGran(ssmid, mymo, start_time, days, coeff)
        if DEBUG:
            print ssmid, start_time, days, coeff, timestep, length, ngran
    else:
        timestep = 0.03125
        length = 2.0
        ngran = 64
    doOneRecursiveSegment(ssmid, mymo, start_time, days, coeff, timestep,
                          length, ngran,  multiplier, CoeffFile, ResidualSumfile,
                          Failedfile,  inputfilename)


def main(argv):
    inputfilepath = argv[0]
    start_time = float(argv[1])
    days = int(argv[2])
    coeff = int(argv[3])
    totaldays = int(argv[4])

    print 'Generating and Fitting Ephems starting:', start_time
    print 'working on file ', inputfilepath
    print 'timespan in days ', days
    print 'number of coefficients ', coeff

    # open output files
    inputfilename = inputfilepath.split("/")
    CoeffFile = open(inputfilename[-1] + '.coef_vartime_' + str(coeff) + '.dat', 'w')
    ResidualSumfile = open(inputfilename[-1] + '.resid_sum_vartime_' + str(coeff) + '.dat', 'w')
    Failedfile = open(inputfilename[-1] + '.failed_' + str(coeff) + '.dat', 'w')

    # get input
    orbit = np.loadtxt(inputfilepath, comments='!!', usecols=(2, 3, 4, 5, 6, 7, 8, 9),
                       delimiter=None, dtype=np.float64, unpack=True)
    ssmid = np.loadtxt(inputfilepath, comments='!!', usecols=(0,), delimiter=None, dtype=np.str, unpack=True)
    q = orbit[0]
    e = orbit[1]
    inc = orbit[2]
    omega = orbit[3]
    argperi = orbit[4]
    t_p = orbit[5]
    H = orbit[6]
    t_0 = orbit[7]

    oo.pyoorb.oorb_init(ephemeris_fname="")
    print 'total days ', totaldays

    # Make Multiplier Dict:
    VMAG_COEFF = 9
    DIST_COEFF = 5
    SE_COEFF = 6

    # Precompute multiplier because
    # we don't want to invert a matrix for every segment
    nPoints = 64
    multipliers = {}
    multipliers['POSITION'] = cg.makeChebMatrix(nPoints + 1, coeff, weight=0.16)
    multipliers['VMAG_X'] = cg.makeChebMatrixOnlyX(nPoints + 1, VMAG_COEFF)
    multipliers['DIST'] = cg.makeChebMatrix(nPoints + 1, DIST_COEFF, weight=0.16)
    multipliers['DIST_X'] = cg.makeChebMatrixOnlyX(nPoints + 1, DIST_COEFF)
    multipliers['SE_X'] = cg.makeChebMatrixOnlyX(nPoints + 1, SE_COEFF)

    # check if only 1 row
    theShape = orbit.shape
    if len(theShape) == 2:
        datalen = theShape[1]
    else:
        datalen = 1
    if DEBUG:
        print 'datalen', datalen
    for i in range(datalen):
        if datalen == 1:
            id = ssmid
            mymo = mo.MovingObject(q, e, inc, omega, argperi, t_p, t_0, objid=ssmid, magHv=H)
        else:
            id = ssmid[i]
            mymo = mo.MovingObject(q[i], e[i], inc[i], omega[i], argperi[i],
                                   t_p[i], t_0[i], objid=ssmid[i], magHv=H[i])

        tmpStartTime = start_time
        while tmpStartTime < start_time + totaldays:
            doOneMonth(id, mymo, tmpStartTime, days, coeff, multipliers, CoeffFile,
                       ResidualSumfile, Failedfile, inputfilename[-1])
            tmpStartTime += days

    CompletedNotice = open(inputfilename[-1] + '.done.txt', 'w')
    print >>CompletedNotice, "Success"


if __name__ == "__main__":
    main(sys.argv[1:])
