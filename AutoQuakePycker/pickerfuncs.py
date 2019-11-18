#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Various functions for the Kurtosis picker of AutoQuakePycker

@author: Stephen Hicks
"""
import numpy as np


def CF_kurtosis(kurt_win_size, tr):
    """
    This is the Kurtosis function of Baillard et al. (2014, SRL)
    (Eq. 10). Sets values outside of window to NaN.

    Parameters
    ----------
    kurt_win_size : kurtosis window length in seconds
    tr : ObsPy trace object

    Returns
    ----------
    cf : kurtosis characteristic function time series
    """
    from scipy.stats import kurtosis

    nwind = int(np.ceil(kurt_win_size * tr.stats.sampling_rate))
    cf = np.zeros(tr.stats.npts)
    for i in range(nwind, tr.stats.npts-1):
        data = tr.data[i-nwind+1: i+1]
        kurt = kurtosis(data, axis=0, bias=False, fisher=False)
        cf[i] = kurt
    cf[0:nwind] = np.nan
    return cf


def kurt_transform_f2(f1, kurt_win_size, tr):
    """
    Apply transformations as per Baillard (BSSA).
    Remove neg slope of the CF from the sliding kurtosis fonction.

    Parameters
    ----------
    f1 : f1 transform time series
    kurt_win_size : kurtosis window length in seconds
    tr : ObsPy trace object

    Returns
    ----------
    f2 : f2 kurtosis transform time series
    """
    nwind = int(np.ceil(kurt_win_size * tr.stats.sampling_rate))
    f2 = np.zeros(len(f1))
    f2[nwind] = f1[nwind]
    # Eq. 11
    for k in range(nwind, len(f1)-1):
        dfi = f1[k+1] - f1[k]
        if dfi >= 0:
            delt = 1
        else:
            delt = 0
        f2[k+1] = f2[k] + delt * dfi
    f2[0:nwind] = np.nan
    f2[-nwind:] = np.nan
    return(f2)


def kurt_transform_f3(f2, kurt_win_size, tr):
    """
    Remove the linear trend of F2.
    Eq. 12

    Parameters
    ----------
    f2 : f2 transform time series
    kurt_win_size : kurtosis window length in seconds
    tr : ObsPy trace object

    Returns
    ----------
    f3 : f3 kurtosis transform time series
    """
    nwind = int(np.ceil(kurt_win_size * tr.stats.sampling_rate))
    a = (f2[len(f2)-nwind-1] - f2[nwind]) / (len(f2)-nwind-1)
    b = f2[nwind]
    f3 = np.zeros(len(f2))
    for k in range(nwind, len(f2)):
        f3[k] = f2[k] - (a*(k-1)+b)
    return f3


def kurt_transform_f4(f3, kurt_win_size, tr):
    """
    Find greatest minima correspond to the greatest onset strengths.

    Parameters
    ----------
    f3 : f3 transform time series
    kurt_win_size : kurtosis window length in seconds
    tr : ObsPy trace object

    Returns
    ----------
    f4 : f4 kurtosis transform time series
    """
    nwind = int(np.ceil(kurt_win_size * tr.stats.sampling_rate))
    f4 = np.zeros(len(f3))
    T = np.zeros(len(f3))
    for i in range(nwind, (len(f3)-1)):
        T[i] = f3[i] - f3[i+1]
        if T[i] < 0:
            f4[i] = T[i]
        elif T[i] >= 0:
            f4[i] = 0.0
        f4[0:nwind] = np.nan
    return(f4)


def get_best_picks(sta, max_wind, orig_time, comps, sptype, cfg):
    """
    Define best Picks
    Choose whether to use fixed min S-P time or variable S-P time depending on
    P-arrival time and Vp/Vs ratio (Wadati method)
    This is quite messy. Could benefit from some improvements.

    Parameters
    ----------
    sta : Station dictionary containing picking results
    max_wind : Window length around theoretical arrival time to select pick in
    orig_time : UTCDateTime

    Returns
    ----------
    f4 : f4 kurtosis transform time series
    """
    best_picks = {}
    best_picks["P"] = [None] * len(cfg.picking.CF_MEAN_SMOOTH_WIND)
    best_picks["S"] = [None] * len(cfg.picking.CF_MEAN_SMOOTH_WIND)
    picks, all_P, all_S, best_p, best_s, = ([] for i in range(5))  # Set empty
    for nphase, phase in enumerate(["P", "S"]):
        pred = sta["pred"][phase]
        for cmp in sta[phase]["picks"]["poss_obs"][phase].keys():
            for nsm, smooth_wind in enumerate(cfg.picking.CF_MEAN_SMOOTH_WIND):
                best_picks[phase][nsm] = []
                if sta[phase]["picks"]["poss_obs"][phase][cmp][nsm]:
                    picks = sta[phase]["picks"]["poss_obs"][phase][cmp][nsm]
                    if phase == "P":
                        picks = [p for p in picks if
                                 np.abs(p[0]-pred) <= max_wind]
                    elif phase == "S" and best_p:
                        if sptype == "const":
                            minsp = cfg.picking.MIN_SP_TIME
                        elif sptype == "dist":
                            minsp = (1.60 - 1) * (best_p[0] - orig_time)
                        picks = [p for p in picks if
                                 np.abs(p[0]-pred) <= max_wind
                                 and p[0] - best_p[0] > minsp]
                if picks:
                    best_loc = picks[0]
                    for pick in sorted(picks, key=lambda x: x[0]):
                        if nsm == 0 and pick[1] <= best_loc[1]:
                            best_loc = pick
                            best_picks[phase][0] = pick
                            lastnsm = nsm
                        elif (nsm > 0 and pick[1] <= best_loc[1]):
                            if (np.abs(pick[0] - best_picks[phase][lastnsm][0])
                                    <= cfg.picking.MAX_DRIFT):
                                best_loc = pick
                                best_picks[phase][nsm] = pick
                                lastnsm = nsm
                    if not best_picks[phase][nsm]:  # Use best from prev smooth
                        best_picks[phase][nsm] = best_picks[phase][nsm-1]
                elif not picks and nsm == 0:  # No pick from smoothest level
                    best_picks[phase][nsm] = []
                    break
                else:
                    best_picks[phase][nsm] = []
        # Now only look at the picks using the least smoothed data.
        if best_picks[phase][-1]:
            if phase == "P":
                all_P.append(best_picks[phase][-1])
                best_p = min(all_P, key=lambda x: x[0])
            elif phase == "S":
                all_S.append(best_picks[phase][-1])
                # If using ZRT system, then prioritise pick on tranverse.
                # Otherwise, take the radial pick.
                if phase == "S" and "R" in comps["S"] and "T" in comps["S"]:
                    all_S_T = [p for p in all_P if p[2][2:3] == "T"]
                    if all_S_T:
                        best_s = min(all_S_T, key=lambda x: x[0])
                    else:
                        best_s = min(all_S, key=lambda x: x[0])
                else:
                    best_s = min(all_S, key=lambda x: x[0])
        else:
            if phase == "P":
                best_p = []
            elif phase == "S":
                best_s = []
    sta["best_obs"] = [best_p, best_s]
    return sta


# def PSswavefilter(rectilinearity, incidence):
#    """
#    p and s wave filter arrays based on Ross et al, 2014 (GJI),  2016
#    input :
#        - rectilinearity  from polarisation analysis
#        - incidence from polarisation analysis (Flinn, 1988)
#    output :
#        - pfilter: p filter (convolution with the raw signal will remove S)
#        - sfilter: s filter (convolution with the raw signal will remove P)
#    """
#    sfilter = (rectilinearity*(np.ones(len(incidence)) -
#                               np.cos(incidence*math.pi/180)))
#    pfilter = rectilinearity*(np.cos(incidence*math.pi/180))
#    return sfilter, pfilter
#
# def polarisation_filters(st_pol, t_new):
#    # 3s window is quoted by Ross et al. (2014, GJI)
#    st_pol.filter("lowpass", freq=3)
#    pol = polarization_analysis(
#            st_pol, 4, 0.01, 4, 6, st_pol[0].stats.starttime+5,
#            st_pol[0].stats.endtime, verbose=False, method='flinn',
#            var_noise=0.0)
#    rectlinearity = pol['rectilinearity']
#    # Shift - Not sure why this is needed to align the real data & pol. filts
#    incidence = pol['incidence']
#    t = pol['timestamp']
#    rectlinearity2 = Trace(data=rectlinearity)
#    rectlinearity2.stats.starttime = t[0]
#    rectlinearity2.stats.delta = t[1]-t[0]
#    t = rectlinearity2.times(reftime=orig_time)
#    r = rectlinearity2.data
#    incidence2 = Trace(data=incidence)
#    incidence2.stats.starttime = t[0]
#    incidence2.stats.delta = t[1]-t[0]
#    i = incidence2.data
#    r_new = np.interp(t_new, t, r)
#    i_new = np.interp(t_new, t, i)
#    alpha = 1
#    DR = np.array([r_new[i]*np.sign(alpha * math.sin((90-i_new[i]))-r_new[i])
#                   for i in
#                  range(len(r_new))])
#    sfilter, pfilter = PSswavefilter(r_new, i_new)
#    return (pfilter, sfilter, DR)
