"""
Compute Chebyshev polynomial coefficients from ephemeris (.ephem) file

Usage:
python generate_coefficients.py <filename> <days> <rows_per_day> <coefficients> <range>

Example:
python generate_coeff_static_tai.py S1.01.xxa.50813.16.ephem.dat 30 32 14 2

Assume that MJDs in .ephem input file are in UTC and convert to TAI.
"""

import numpy as np
import chebfit as cg
import sys


def utc2tai(tUTC):
    if tUTC < 49169:
        raise RuntimeError('Function not defined for UTC < 49169')
    elif tUTC < 49534:
        return tUTC + 0.000324074074  # 28s
    elif tUTC < 50083:
        return tUTC + 0.000335648148  # 29s
    elif tUTC < 50630:
        return tUTC + 0.000347222222  # 30s
    elif tUTC < 51179:
        return tUTC + 0.000358796296  # 31s
    elif tUTC < 53736:
        return tUTC + 0.00037037037   # 32s
    elif tUTC < 54832:
        return tUTC + 0.000381944444  # 33s
    else:
        return tUTC + 0.000393518518  # 34s


vUtc2Tai = np.vectorize(utc2tai)


def get_coeffs_position(t, ra, dec, draskydt, ddecdt, ng, npo, coeff):
    p, dec_resid, dec_rms = cg.chebfit(t, dec, ddecdt, nPoly=coeff)
    rap, ra_coord_resid, ra_rms = cg.chebfit(t, ra, draskydt/np.cos(np.pi*dec/180.), nPoly=coeff)
    ra_sky_resid = (ra_coord_resid)*np.cos(np.pi*dec/180.)
    return p, rap, 3600.0*1000.0*np.max(np.sqrt(dec_resid**2 + ra_sky_resid**2))


def get_coeffs_vmag(t, vmag, ng, npo, coeff):
    p, resid, rms = cg.chebfit(t, vmag, None, nPoly=coeff)
    return p, np.max(np.abs(resid))


def get_coeffs_se(t, se, ng, npo, coeff):
    p, resid, rms = cg.chebfit(t, se, None, nPoly=coeff)
    return p, np.max(np.abs(resid))


def get_coeffs_dist(t, dist, distdt, ng, npo, coeff):
    p, resid, rms = cg.chebfit(t, dist, distdt, nPoly=coeff)
    return p, np.max(np.abs(resid))


def three_sixy_to_neg(element, min, max):
    if (min < 10) & (max > 350):
        if element > 180:
            return element - 360.
        else:
            return element
    else:
        return element


def main(argv):
    inputfilepath = argv[0]
    days_per_object = argv[1]
    rows_per_day = int(argv[2])
    coeff = int(argv[3])
    range = np.float(argv[4])

    short = False
    print 'working on file ', inputfilepath
    print 'timespan in days ', days_per_object
    print 'rows per day ', rows_per_day
    print 'number of coefficients ', coeff

    # SETTINGS #
    daystart = int(range*rows_per_day)
    skip = 1

    ###################################
    # open output files
    inputfilename = inputfilepath.split("/")
    CoeffFile = open(inputfilename[-1] + '.coef_vartime_' + str(coeff) + '.dat', 'w')
    ResidualSumfile = open(inputfilename[-1] + '.resid_sum_vartime_' + str(coeff) + '.dat', 'w')

    # utc2tai = 32./86400.
    vUtc2Tai = np.vectorize(utc2tai)
    # get input
    ephem = np.loadtxt(inputfilepath, comments='#', usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 14, 15),
                       delimiter=None, dtype=np.float64, unpack=True)
    if short:
        ssmid = ephem[0]
        t = vUtc2Tai(ephem[1])
        ra = ephem[3]
        dec = ephem[4]
        draskydt = ephem[6]
        ddecdt = ephem[7]
        dist = ephem[2]
        distdt = ephem[5]
        vmag = ephem[8]
        se = ephem[10]
    else:
        ssmid = ephem[0]
        t = vUtc2Tai(ephem[2])
        ra = ephem[4]
        dec = ephem[5]
        draskydt = ephem[7]
        ddecdt = ephem[8]
        dist = ephem[3]
        distdt = ephem[6]
        vmag = ephem[9]
        se = ephem[11]

    advance_next_day = int(rows_per_day)*int(days_per_object) + 1
    objectcount = 1
    outerday0 = 0
    rowarray = np.array([])
    residarray = np.array([])
    deltaresidarray = np.array([])
    vmagresidarray = np.array([])
    v_360_to_neg = np.vectorize(three_sixy_to_neg)

    while outerday0 < len(ephem[0]):
        day0 = outerday0
        day1 = outerday0 + daystart
        rows = 1
        # NEW OBJECT
        # For one object over the course of 1 month or 43200 minutes
        while day0 < advance_next_day*objectcount - 1:
            ngran = day1 - day0
            npoint = day1 - day0
            pdec, pra, p_resid = get_coeffs_position(t[day0:day1 + 1:skip],
                                                     v_360_to_neg(ra[day0:day1 + 1:skip],
                                                                  np.min(ra[day0:day1 + 1:skip]),
                                                                  np.max(ra[day0:day1 + 1:skip])),
                                                     dec[day0:day1 + 1:skip], draskydt[day0:day1 + 1:skip],
                                                     ddecdt[day0:day1 + 1:skip], ngran, npoint, coeff)
            d, d_resid = get_coeffs_dist(t[day0:day1 + 1:skip], dist[day0:day1 + 1:skip],
                                         distdt[day0:day1 + 1:skip], ngran, npoint, 5)
            v, v_resid = get_coeffs_vmag(t[day0:day1 + 1:skip], vmag[day0:day1 + 1:skip], ngran, npoint, 9)
            s, s_resid = get_coeffs_se(t[day0:day1 + 1:skip], se[day0:day1 + 1:skip], ngran, npoint, 8)
            vmagresidarray = np.append(vmagresidarray,  v_resid)
            deltaresidarray = np.append(deltaresidarray,  d_resid)
            residarray = np.append(residarray, p_resid)


            print >>ResidualSumfile, "%i %i %.14f %.14f %.14f %.14e %.14e %.14e %.14e %s"%(ssmid[day0], rows, t[day0], t[day1], t[day1] - t[day0], p_resid, d_resid, v_resid, s_resid, inputfilename[-1])
            print >>CoeffFile, "%i %i %.6f %.6f %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e %.14e"%(0,ssmid[day0], t[day0], t[day1], pra[0],pra[1],pra[2],pra[3],pra[4],pra[5],pra[6],pra[7],pra[8],pra[9],pra[10],pra[11],pra[12],pra[13],pdec[0],pdec[1],pdec[2],pdec[3],pdec[4],pdec[5],pdec[6],pdec[7],pdec[8],pdec[9],pdec[10],pdec[11],pdec[12],pdec[13],d[0],d[1],d[2],d[3],d[4],v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],s[0],s[1],s[2],s[3],s[4],s[5])

            # advance to the next day if less than 6 points left in month
            day0 = day1
            day1 = day0 + daystart
            rows = rows + 1
        rowarray = np.append(rowarray, rows)
        objectcount = objectcount + 1
        outerday0 = outerday0 + advance_next_day

    print np.min(rowarray), np.max(rowarray), np.mean(rowarray)
    print np.min(residarray), np.max(residarray), np.mean(residarray)
    print np.min(deltaresidarray), np.max(deltaresidarray), np.mean(deltaresidarray)
    print np.min(vmagresidarray), np.max(vmagresidarray), np.mean(vmagresidarray)
    CompletedNotice = open(inputfilename[-1] + '.done.txt', 'w')
    print >>CompletedNotice, "Success"

if __name__ == "__main__":
    main(sys.argv[1:])
