#!/usr/bin/env python
"""
main_script.py
====================================
test
30/09/19
Tried using polarisation analysis as well (e.g. Ross et al. and Bailard et al.)
- P and S filters and Dip-Rectilinearity value but did not produce consistent
results - bad polarisations - bad rotations - noisy data???
These sections have been commented out.
04/11/19
Found that separate filters are needed for P and S-waves
Added integration for S-wave picking

Final refinement takes S-picks on the T-component only to ensure that S is not
mispicked as Sp precursors.
Also computes magnitudes.
TODO: Consider picking on radial as well as transverse? (Preferentially select tranverse)
TODO: Check memory usage - are some variables being appended to?
TODO: Process TRN events
TODO: 8/11 Discovered some duplicates with Lidong's catalog (auto origin time
from lassie was too far off to be detected
TODO: Implement Sippl approach: If RMS too high, progressively remove picks with largest 
residuals until suitable RMS is found
TODO: On first step, only use stations that are closest, then incorporate more distant 
arrivals in later iteration. Do this 1 pick at a time if after relocation the single residual
of the new pick fell below a distance‐dependent goal residual (P:0.6 s + (traveltime)*0.015);
(S: 1.0 s + (traveltime)*0.015
"""

from obspy import read_events, read, Stream, read_inventory
from obspy.core.event import (
        Pick, WaveformStreamID, EventDescription, QuantityError, Magnitude,
        Amplitude, Catalog)
from obspy.core.event.header import EvaluationMode
from obspy.geodetics.base import gps2dist_azimuth
from obspy.signal.trigger import classic_sta_lta
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import numpy as np
from scipy.stats import kurtosis, trim_mean
import os
import subprocess
import glob
import time
import multiprocessing
import more_itertools as mit
import shutil
import sys
import logging
from itertools import product

logger = multiprocessing.log_to_stderr(logging.DEBUG)



# Start of parameters to define

# INPUT OPTIONS
lassie_cat_file = "../new_event_detections.xml"
DIR_TO_EVENTDIRS = "../events/events"

# OUTPUT OPTIONS
FORCE_RECALC = False  # Choose whether or not to overwrite existing events

# SNR / NOISE OPTIONS
SNR_WIN = [5, 1]  # Long-time short time window for finding trace noise level
SNR_FREQ = [4, 8]  # Bandpass [low-cut, high-cut] filter freq for computing SNR
MAX_NOISE_LEVEL = 50.0  # Max. noise level to reject trace
SNR_SMOOTH_WIN = 1.0  # Window length for SNR smoothing

# STATION SELECTION OPTIONS
MAX_DIST = 275  # Max station-event distance to pick on station
STA_BLACKLIST = ["BAUV", "TRNT"]

# KURTOSIS PICKING OPTIONS
# Components to use at different refinement stages:
CMPS_REFINE_1 = {"P": ["Z", "H"], "S": ["1", "2", "N", "E"]}
CMPS_REFINE_2 = {"P": ["Z", "H"], "S": ["T", "R"]}
INTEGRATE_S = True  # Choose whether to integrate to pick S-wave arrival
# Maximum allowed difference between theoretical arrival time & pick (in secs):
# (smaller values = faster processing)
# (final value should be set to maximum possible travel-time residual)
MAX_PICK_DIFF_REFINE1 = 12.0
MAX_PICK_DIFF_REFINE2 = 3.5
KURT_WINS = [2, 4, 6]  # Kurtosis smoothing window lengths (s) for f1
# P[0]- & S[1]-wave (no. filters need to be the same):
FILT_WINS = {"P": [[4, 12], [10, 18]], "S": [[2, 6], [3, 8]]}
CF_MEAN_SMOOTH_WIND = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.08, 0.06]
MAX_DRIFT = 0.12  # Max dt between increasingly roughened windows of f4
MIN_SP_TIME = 2.5  # Minimum S-P time in seconds for first refinement
# Most-smoothed Kurtosis values (f1) to choose pick weighting [0, 1, 2, 3, 4]
KURT2WGHT = {"P": [6, 4.8, 4, 2.9, 0], "S": [6, 4.8, 4, 2.9, 0]}
# Translation between weight and arrival time uncertainty for NonLinLoc
# [0, 1, 2, 3, 4]
T_ERRORS = [0.1, 0.2, 0.6, 1.2, 99.9]

# MAGNITUDE CALCULATION OPTIONS
ML_CAL = [1., 3.01e-3, 3.]  # ML calibration abc values (default: BakunJoyner)
MAX_DIST_ML = 400.0  # Max epicentral distance to compute station ML value
V_MIN = 3.0  # Minimum velocity to compute amplitude pick window
MAX_WINDOW = 150.0  # Maximum window legnth for amplitude picking
MAG_HIGHPASS = 2  # High-pass filter for amplitude picking

# PLOTTING OPTIONS
DO_PLOT_1 = True
DO_PLOT_2 = True
PLOT_WIND_PICK_DIFF = 4.0  # Secs to plot around picked arr.
# Larger value good for initial debugging / tuning parameters

# End of parameters to define


def add_picks(sta, e, orig_time, st):
    """
    Add picks made to ObsPy event structure.
    Args:
        sta (dict): variable station & picking data.
        e (ObsPy event object): event.
        orig_time (ObsPy UTCDateTime object): event origin time.
        st (Obspy stream object): stream of ObsPy traces containing data.
    Returns:
        e (ObsPy event object): event containg added picks.

    """
    for n_p, p in enumerate(sta["best_obs"]):
        if n_p == 0:
            phase = "P"
        elif n_p == 1:
            phase = "S"
        if p and p[3] != 4:
            pick = Pick(waveform_id=WaveformStreamID(
                network_code=st[0].stats.network,
                station_code=st[0].stats.station,
                channel_code=p[2]),
                time=p[0], phase_hint=phase,
                evaluation_mode=EvaluationMode("automatic"),
                time_errors=QuantityError(uncertainty=T_ERRORS[p[3]]))
            e.picks.append(pick)
    return(e)


def compute_magnitude(evt, st):
    """
    Compute local magnitude for event. Uses the same approach as Bie et al.
    (2019, SRL), although we don't bother with checking signal-to-noise ratio
    this time.
    Args:
        
    """
    sta_amps = []
    new_origin = evt.preferred_origin()
    stations = [pick.waveform_id.station_code for pick in evt.picks if
                pick.phase_hint == "P"]
    for station in stations:
        PAZ_WA = {'poles': [-6.283 + 4.7124j, -6.283 - 4.7124j],
                  'zeros': [0 + 0j], 'gain': 1.0, 'sensitivity': 2080}
        st_sta = st.select(station=station)
        st_mag = st_sta.copy()
        for tr in st_mag.select(component="H"):
            st_mag.remove(tr)
        for tr in st_mag:
            if len(tr.data) == 0:
                st_mag.remove(tr)
        inv = read_inventory("../../metadata/{:}.{:}.xml".format(
                st_mag[0].stats.network, station)).select(station=station)
        st_mag.detrend("demean")
        st_mag.detrend("linear")
        st_mag.taper(max_percentage=0.05)
        try:
            st_mag.remove_response(inv, output="VEL",
                                   pre_filt=[0.05, 0.08, 18, 19])
        except ValueError as e:
            print(e, station, st_mag)
        st_mag.simulate(paz_remove=None, paz_simulate=PAZ_WA, water_level=60.0)
        st_mag.filter("highpass", freq=MAG_HIGHPASS)
        pick = [pick for pick in evt.picks if
                pick.waveform_id.station_code == station
                and pick.phase_hint == "P"]
        if pick:
            pick = pick[0]
        arrival = [arrival for arrival in new_origin.arrivals
                   if arrival.pick_id == pick.resource_id][0]
        dist_km = arrival.distance * 111
        depth = new_origin.depth / 1000
        hyp_dist = np.sqrt(dist_km**2 + depth**2)
        if dist_km >= MAX_DIST_ML:
            continue
        Lg_end = new_origin.time + dist_km / V_MIN + 30
        if Lg_end - new_origin.time >= MAX_WINDOW:
            Lg_end = new_origin.time + MAX_WINDOW
        st_mag.trim(pick.time-0.25, Lg_end)
        st_amps = []
        st_amps_t = []
        for tr in st_mag:
            st_amps.append(np.max(np.abs(tr.data)))
            st_amps_t.append(tr.times(reftime=new_origin.time)
                             [np.argmax(np.abs(tr.data))])
        amp_sta = np.max(st_amps)
        n_max = np.argmax(st_amps)

        amp = Amplitude(generic_amplitude=amp_sta, type="AML", unit="m",
                        waveform_id=WaveformStreamID(
                                network_code=st_mag[n_max].stats.network,
                                station_code=st_mag[n_max].stats.station,
                                location_code=st_mag[n_max].stats.location,
                                channel_code=st_mag[n_max].stats.channel),
                        magnitude_hint="ML", pick_id=arrival.pick_id
                        )
        amp = amp_sta * 1000
        ML = (
                np.log10(amp) + ML_CAL[0] * np.log10(hyp_dist/100) +
                ML_CAL[1]*(hyp_dist-100) + ML_CAL[2]
                )
        sta_amps.append(ML)
    ML_ave = trim_mean(sta_amps, 0.25)
    evt.magnitudes.append(Magnitude(
            mag=ML_ave, magnitude_type="ML_voila",
            origin_id=new_origin.resource_id,
            station_count=len(sta_amps)))
    evt.preferred_magnitude_id = evt.magnitudes[0].resource_id
    evt.preferred_magnitude_id = evt.magnitudes[0].resource_id
    return(evt)


def compute_noise_levels(tr):
    """
    Get simple estimate of trace noise level using STA/LTA function - similar
    to thhat described by Baillard et al. (2014, BSSA).
    """
    tr_snr = tr.copy()
    tr_snr.filter("bandpass", freqmin=SNR_FREQ[0], freqmax=SNR_FREQ[1])
    wa = int(SNR_WIN[1]*tr.stats.sampling_rate)
    wb = int(SNR_WIN[0]*tr.stats.sampling_rate)
    if len(tr_snr.data) < wa or len(tr_snr.data) < wb:  # Prevent failing due to en(data) < nlta error
        noise_level = 100.0
        return noise_level
    snr = classic_sta_lta(tr_snr.data, wa, wb)
    snr_smooth = do_smooth(snr, SNR_SMOOTH_WIN, tr.stats.sampling_rate)
    thresh_snr = np.nanmax(snr_smooth) * 0.4
    A = (snr_smooth - thresh_snr)
    A = A[np.where(A > 0)]
    if len(snr_smooth[wb:-wa]) == 0:  # In case zerodivision error
        noise_level = 9999.9
        return noise_level
    noise_level = (len(A) / len(snr_smooth[wb:-wa])) * 100
    return noise_level


def copytree(src, dst, symlinks=False, ignore=None):
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            shutil.copytree(s, d, symlinks, ignore)
        else:
            shutil.copy2(s, d)


def CF_kurtosis(kurt_win_size, tr):
    """
    This is the Kurtosis function of Baillard et al. (2014, SRL)
    (Eq. 10). Sets values outside of window to NaN.
    """
    nwind = int(np.ceil(kurt_win_size * tr.stats.sampling_rate))
    cf = np.zeros(tr.stats.npts)
    for i in range(nwind, tr.stats.npts-1):
        data = tr.data[i-nwind+1: i+1]
        kurt = kurtosis(data, axis=0, bias=False, fisher=False)
        cf[i] = kurt
    cf[0:nwind] = np.nan
    return cf


def get_best_picks(sta, max_wind, orig_time, comps, sptype):
    """
    Define best Picks
    Choose whether to use fixed min S-P time or variable S-P time depending on
    P-arrival time and Vp/Vs ratio (Wadati method)
    This is quite messy. Could benefit from some improvements.
    """
    best_picks = {}
    best_picks["P"] = [None] * len(CF_MEAN_SMOOTH_WIND)
    best_picks["S"] = [None] * len(CF_MEAN_SMOOTH_WIND)
    picks, all_P, all_S, best_p, best_s, = ([] for i in range(5))  # Set empty
    for nphase, phase in enumerate(["P", "S"]):
        pred = sta["pred"][phase]
        for cmp in sta[phase]["picks"]["poss_obs"][phase].keys():
            for nsm, smooth_wind in enumerate(CF_MEAN_SMOOTH_WIND):
                best_picks[phase][nsm] = []
                if sta[phase]["picks"]["poss_obs"][phase][cmp][nsm]:
                    picks = sta[phase]["picks"]["poss_obs"][phase][cmp][nsm]
                    if phase == "P":
                        picks = [p for p in picks if
                                 np.abs(p[0]-pred) <= max_wind]
                    elif phase == "S" and best_p:
                        if sptype == "const":
                            minsp = MIN_SP_TIME
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
                                    <= MAX_DRIFT):
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


def get_theoretical_tt(event):
    """
    Find theoretical arrival times from NLLoc Time2EQ (containing 1 event).
    """
    e = event.copy()
    orig = e.origins[0]
    ev_id = e.event_descriptions[0].text
    e.picks = []
    eqsrce = "EQSRCE {:} LATLON {:7.4f} {:7.4f} {:6.2f} 0.0\n".format(
            ev_id, orig.latitude, orig.longitude,  orig.depth/1000)
    f = open("run.in", "r")
    w = open("run.tmpsynth", "w")
    for l in f:
        w.write(l)
    w.write(eqsrce)
    w.close()
    p1 = subprocess.Popen("./1_compute_theor_arrtimes_NLLOC.sh", shell=True)
    p1.communicate()
    ff = open("out_pred_tt.dat", "r")
    for n, l in enumerate(ff):
        if n > 0 and l[0:1] != "#":
            phase = l.split()[4]
            arr_time = (orig.time
                        + float(l.split()[7][-2:].lstrip()) * 60
                        + float(l.split()[8]))
            sta_code = l.split()[0]
            e.picks.append(Pick(
                    time=arr_time, phase_hint=phase,
                    waveform_id=WaveformStreamID(station_code=sta_code)))
    ff.close()
    return e.picks


def kurt_transform_f2(f1, kurt_win_size, tr):
    """
    Apply transformations as per Baillard (BSSA).
    Remove neg slope of the CF from the sliding kurtosis fonction.
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
    """
    nwind = int(np.ceil(kurt_win_size * tr.stats.sampling_rate))
    a = (f2[len(f2)-nwind-1] - f2[nwind]) / (len(f2)-nwind-1)
    b = f2[nwind]
    f3 = np.zeros(len(f2))
    for k in range(nwind, len(f2)):
        f3[k] = f2[k] - (a*(k-1)+b)
    return f3


def process_events(cat_data, n_run):
    """
    This is the main function that runs a given set of events and does the
    multi-stage refinement.
    """
    if FORCE_RECALC is True:
        w = open("refined_events.dat", "w")
        w.close()
    # Prepare directory
    if (os.path.exists("runs/run{:}".format(n_run))
            and os.path.isdir("runs/run{:}".format(n_run))):
        shutil.rmtree("runs/run{:}".format(n_run))
    copytree("NLLOC_run", "runs/run{:}".format(n_run))
    os.chdir("runs/run{:}".format(n_run))
    for n_ev, ev in enumerate(cat_data):
        start = time.time()
        ev_id = ev.event_descriptions[0].text
        sys.stdout.flush()
        ev_dict = {}
        ev_dict["stations"] = {}
        orig_lat, orig_lon = [ev.origins[0].latitude, ev.origins[0].longitude]
        logging.debug("startint logging")
        st = read("../../{:}/{:}/MSEED/*.msd".format(DIR_TO_EVENTDIRS, ev_id,
                   format="MSEED"))
            
        print(n_run, ev_id)
        for n_tr, tr in enumerate(st):
            if st[n_tr].stats.sampling_rate > 40.0:
                try:
                    st[n_tr].resample(40)
                except ZeroDivisionError:
                    continue
        st1, st2, st_mag = [st.copy(), st.copy(), st.copy()]
        # Append distance to trace
        stations_data = sorted(set([tr.stats.station for tr in st if not
                                    tr.stats.station in STA_BLACKLIST]))
        stations_dist = {sta_code: gps2dist_azimuth(
            sta_locs[sta_code]["lat"], sta_locs[sta_code]["lon"],
            orig_lat, orig_lon)[0] for sta_code in stations_data
            if gps2dist_azimuth(
                sta_locs[sta_code]["lat"], sta_locs[sta_code]["lon"],
                orig_lat, orig_lon)[0]/1000 <= MAX_DIST}
        path_to_figs = "../../{:}/{:}/figs".format(DIR_TO_EVENTDIRS, ev_id)
        if not os.path.exists(path_to_figs):
            os.mkdir(path_to_figs)
        print("Doing first refinement")
        sys.stdout.flush()
        if "R" in CMPS_REFINE_1["S"] or "T" in CMPS_REFINE_1["S"]:
            rot = True
        else:
            rot = False
        evt_refine_1, rms, found = refine_events(
                st1, stations_dist, CMPS_REFINE_1, MAX_PICK_DIFF_REFINE1, ev,
                DO_PLOT_1, 1, fig, "const", path_to_figs, ev_dict, ev_id, rot
                )
        if found is False:
            continue
        print("RMS = ", rms)
        sys.stdout.flush()
        prev_rms = rms
        print("Doing second refinement")
        sys.stdout.flush()
        if "R" in CMPS_REFINE_2["S"] or "T" in CMPS_REFINE_2["S"]:
            rot = True
        else:
            rot = False
        evt_refine_2, rms, found = refine_events(
                st2, stations_dist, CMPS_REFINE_2, MAX_PICK_DIFF_REFINE2,
                evt_refine_1, DO_PLOT_2, 2, fig, "dist", path_to_figs, ev_dict,
                ev_id, rot
                )
        if found is False:
            continue
        print("RMS = ", rms)
        if rms > prev_rms * 1.25:
            print("RMS is significantly increasing (*25%) - skipping event")
            continue
        prev_rms = rms
        evt_refine_2 = compute_magnitude(evt_refine_2, st_mag)
        write_evt(evt_refine_2, ev_id)
        end = time.time()
        print("Time taken for event: {:3.1f} mins".format((end-start)/60))


def process_station(sta, st_sta, dist, ev_dict, evt, orig_time, phase, cmps):
    """
    Get ready to process data. Sets up emtpy arrays.
    Checks noise level of data.
    """
    sta["ncha"] = len(st_sta)
    sta["lenD"] = int(st_sta[0].stats.npts)
    sta["samplerate"] = st_sta[0].stats.sampling_rate
    sta["nwinsamp"] = int(np.ceil(KURT_WINS[0] * sta["samplerate"]))
    st_sta.detrend()
    st_sta.detrend(type="demean")
    for n_tr, _ in enumerate(st_sta):
        st_sta[n_tr].stats.distance = dist
    # Set empty arrays for filling during process
    sta["filtwins_check"] = {}
    sta["picks"] = {}
    sta["picks"]["poss_obs"] = {}
    sta["picks"]["poss_obs"]["P"] = {}
    sta["picks"]["poss_obs"]["S"] = {}
    sta["noise_levels"] = np.zeros(len(st_sta))
    for n_tr, tr in enumerate(st_sta):
        sta["noise_levels"][n_tr] = compute_noise_levels(tr)
        if sta["noise_levels"][n_tr] < MAX_NOISE_LEVEL:
            sta = process_trace(n_tr, tr, sta, orig_time, cmps)
    return(sta)


def process_trace(n_tr, tr, sta, orig_time, cmps):
    """
    Process picking on each trace parsed by process_station.
    """
    cmp = tr.stats.channel[2:3]
    sta[cmp] = {}
    sta[cmp]["times"] = tr.times(reftime=orig_time)

    sta[cmp]["tr_results"] = np.zeros(
            (len(FILT_WINS["P"]), sta["lenD"])
            )
    sta[cmp]["f1_results"] = np.zeros(
            (len(FILT_WINS["P"]), len(KURT_WINS), sta["lenD"])
            )
    sta[cmp]["f1_mean"] = np.zeros(sta["lenD"])
    sta[cmp]["f3_results"] = np.zeros(
            (len(FILT_WINS["P"]), len(KURT_WINS), sta["lenD"])
            )
    sta[cmp]["f3_mean_smooth"] = np.zeros(
            (len(CF_MEAN_SMOOTH_WIND), sta["lenD"])
            )
    sta[cmp]["f4_all"] = np.zeros((len(CF_MEAN_SMOOTH_WIND), sta["lenD"]))
    sta[cmp]["f1_mean_smooth"] = np.zeros(sta["lenD"])
    # Get suitable filters (exclude those fully outside Nyquist freq.)
    for phase in ["P", "S"]:
        if cmp in cmps[phase]:
            sta["picks"]["poss_obs"][phase][cmp] = {}
            sta[cmp]["filtwins_check"] = [
                    filt_win for filt_win in FILT_WINS[phase]
                    if filt_win[0] < sta["samplerate"] / 2
                    ]
        if INTEGRATE_S is True:
            tr.integrate()

    for n_filt, filt in enumerate(sta[cmp]["filtwins_check"]):
        # Ensure that filter covers sample rate / 2
        if (tr.stats.sampling_rate / 2) <= filt[0]:
            print("Skipping this Kurtosis run due to sample rate/2<f")
            continue
        tr.filter("bandpass", freqmin=filt[0], freqmax=filt[1])
        try:
            sta[cmp]["tr_results"][n_filt] = tr.data
        except ValueError:  # If input array length is inconsistent
            continue
        # Loop over kurtosis windows
        for n_kurt, kurt_win_s in enumerate(KURT_WINS):
            f1 = CF_kurtosis(kurt_win_s, tr)
            sta[cmp]["f1_results"][n_filt, n_kurt] = f1  # Needed for weights
            f2 = kurt_transform_f2(f1, kurt_win_s, tr)
            f3 = kurt_transform_f3(f2, kurt_win_s, tr)

            sta[cmp]["f3_results"][n_filt, n_kurt] = f3
    sta[cmp]["f1_mean"] = np.nanmean(sta[cmp]["f1_results"], axis=0)[0]
    sta[cmp]["f1_mean_smooth"] = do_smooth(
            sta[cmp]["f1_mean"], CF_MEAN_SMOOTH_WIND[0], tr.stats.sampling_rate
            )
    # ^ Throws up a warning first time due to NaN slices
    # Compute mean CF and final kurtosis transform
    f3_mean = np.nanmean(sta[cmp]["f3_results"], axis=0)[0]

    for nsm, smooth_wind in enumerate(CF_MEAN_SMOOTH_WIND):
        sta[cmp]["f3_mean_smooth"][nsm] = do_smooth(
                f3_mean, smooth_wind, tr.stats.sampling_rate
                )
        f4 = kurt_transform_f4(sta[cmp]["f3_mean_smooth"][nsm],
                               np.max(KURT_WINS), tr)
        sta[cmp]["f4_all"][nsm] = f4

        # Now pick (avoiding end and beginning of signal)
        # Pick the P-waves
        if cmp in cmps["P"]:
            sta["picks"]["poss_obs"]["P"][cmp][nsm] = []
            # Find points where Kurt<0 & doesn't look like S-wave
            p_cands = np.argwhere((f4 < 0.0))
            for idx in p_cands.tolist():
                kurt_wgt = np.min(np.where(np.array(
                        KURT2WGHT["P"] <= sta[cmp]["f1_mean_smooth"][idx])))
                sta["picks"]["poss_obs"]["P"][cmp][nsm].append([
                        orig_time+sta[cmp]["times"][idx][0], f4[idx][0],
                        tr.stats.channel, kurt_wgt, idx,
                        sta[cmp]["times"][idx][0]
                        ])
        # Pick the S-waves
        if cmp in cmps["S"]:
            sta["picks"]["poss_obs"]["S"][cmp][nsm] = []

            # Find points where Kurt<0 & doesn't look like S-wave
            s_cands = np.argwhere((f4 < 0.0))
            for idx in s_cands.tolist():
                kurt_wgt = np.min(np.where(np.array(KURT2WGHT["S"]
                                           <= sta[cmp]["f1_mean_smooth"][idx]))
                                  )
                sta["picks"]["poss_obs"]["S"][cmp][nsm].append([
                        orig_time+sta[cmp]["times"][idx][0], f4[idx][0],
                        tr.stats.channel, kurt_wgt, idx,
                        sta[cmp]["times"][idx][0]
                        ])
    return(sta)


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


def refine_events(st, stations_dist, cmps, max_pick_diff, evt, do_plot,
                  n_refine, fig, sptype, path_to_figs, ev_dict, ev_id,
                  rot=False):
    """
    This is one of the main functions of this code. It processes a refinement
    attempt for each event.
    """
    if do_plot is True:
        pdf = matplotlib.backends.backend_pdf.PdfPages(
                "{:}/refined{:}.pdf".format(path_to_figs, n_refine))
    orig_time = evt.origins[0].time
    ev_dict["picks"] = get_theoretical_tt(evt)
    evt.picks = []   # Empty picks from previous refinement run
    for station, dist in sorted(stations_dist.items(), key=lambda x: x[1]):
        st_sta = st.select(station=station)
        if len(st_sta) >= 3:
            if rot is True:
                st_sta, rot_done = rotate(st_sta, evt)
                if rot_done is False:
                    continue
            sta = {}
            sta["pred"] = {}
            sta["best_obs"] = {}
            for phase in ["P", "S"]:
                sta["pred"][phase] = [
                        p.time for p in ev_dict["picks"]
                        if p.waveform_id.station_code == station and
                        p.phase_hint == phase][0]
                st_proc = Stream()
                # Ensure that trace is padded by longest Kurtosis window
                st_ = st_sta.slice(
                        sta["pred"][phase]-max_pick_diff-np.max(KURT_WINS)-1,
                        sta["pred"][phase]+max_pick_diff+np.max(KURT_WINS))
                for tr in st_:
                    if tr.stats.channel[2:3] in cmps[phase]:
                        st_proc.append(tr)
                if len(st_proc) > 0:
                    sta[phase] = {}
                    sta[phase] = process_station(
                            sta[phase], st_proc, dist, ev_dict,
                            evt, orig_time, phase, cmps)
            try:
                sta = get_best_picks(sta, max_pick_diff, orig_time,
                                     cmps, sptype)
            except:
                continue
            if sta["best_obs"][0]:
                evt = add_picks(sta, evt, orig_time, st_)
                if do_plot is True:
                    plot_sta_results(sta, st_sta, "refine{:}".format(
                            n_refine), orig_time, pdf, fig, cmps)

    if do_plot is True:
        pdf.close()
    evt_refine_out = relocate_nlloc(evt)
    if evt_refine_out == []:
        print("Relocation failed - skipping event")
        found = False
        rms = 100.0
    else:
        evt_refine_out.event_descriptions = [(EventDescription(text=ev_id))]
        rms = evt_refine_out.preferred_origin().quality.standard_error
        found = True
    return(evt_refine_out, rms, found)


def relocate_nlloc(evt_refined):
    """
    Carry out relocation using newly picked events in NonLinLoc.
    """
    evt_refined.write("tmp.nlloc", format="NLLOC_OBS")
    # Now run the location
    os.system("./3_compute_locs.sh")
    if os.path.exists("loc/last.hyp"):
        ev_new = read_events("loc/last.hyp", format="NLLOC_HYP")[0]
    else:
        print("No output NLLOC file for initial relocation\n"
              "Maybe not enough arrivals or something went wrong\n")
        ev_new = []
    return ev_new


def rotate(st, evt):
    """
    Rotate raw 3-component data to ZRT coordinate system.
    Throws up flag if rotation cannot be done.
    """
    done = True
    origin = evt.origins[0]
    inv = (read_inventory("../../metadata/{:}.{:}.xml".format(
            st[0].stats.network, st[0].stats.station)).
           select(station=st[0].stats.station))
    dist, az, baz = gps2dist_azimuth(
            origin.latitude, origin.longitude,
            inv[0][0].latitude, inv[0][0].longitude)
    st_sta_rot = st.copy()
    try:
        st_sta_rot.rotate("->ZNE", inventory=inv).rotate(
            'NE->RT', back_azimuth=baz)
    except ValueError:
        print("For rotation, components don't have the same timespan - "
              "skipping station")
        done = False
    return(st_sta_rot, done)


def do_smooth(d, WT, sample_rate):
    """
    Simple windowing smoothing function.
    (takes data to right of smoothing point).
    """
    d_smooth = np.zeros(len(d))
    Wt = int(np.ceil(sample_rate*WT))
    for i in range(len(d)-Wt):
        d_smooth[i] = np.mean(d[i: i+Wt])
    d_smooth[0:Wt+100] = np.nan  # +100 removes "edge effects" at start of f4
    return(d_smooth)


def kurt_transform_f4(f3, kurt_win_size, tr):
    """
    Find greatest minima correspond to the greatest onset strengths.
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


def plot_record_section(st_sta):
    """
    TODO
    """
    fig = plt.figure()
    st_sta.plot(type="section", orientation="horizontal",
                vertical_scaling_range=5e3, fig=fig)
#    for tr in st_sta:
#    plt.scatter([30], [60])


def plot_sta_results(sta, st_sta, runtype, orig_time, pdf, fig, cmps):
    """
    Plot picking results for each station
    """
    co = ["blue", "orange"]
    ct = ["green", "red"]
    for n_tr, tr in enumerate(st_sta):
        cmp = tr.stats.channel[2:3]
        if cmp in cmps["P"]:
            phase = "P"
            nphs = 0
        elif cmp in cmps["S"]:
            phase = "S"
            nphs = 1
        sta_ = sta[phase]
        try:
            sta_[cmp]
        except:
            continue
        ax1 = plt.subplot(4, len(st_sta), n_tr+1)
        ax2 = plt.subplot(4, len(st_sta), n_tr+len(st_sta)+1, sharex=ax1)
        ax3 = plt.subplot(4, len(st_sta), n_tr+2*(len(st_sta))+1, sharex=ax1)
        ax4 = plt.subplot(4, len(st_sta), n_tr+3*(len(st_sta))+1, sharex=ax1)
        if len(sta_[cmp]["times"]) != len(sta_[cmp]["f1_mean_smooth"]):
            continue
        ax2.plot(sta_[cmp]["times"], sta_[cmp]["f1_mean_smooth"],
                 label="F$_4$")
        for n_filt, filt in enumerate(sta_[cmp]["filtwins_check"]):
            ax1 = ax1.twinx()
            ax1.plot(sta_[cmp]["times"], sta_[cmp]["tr_results"][n_filt],
                     linewidth=0.2, color="black")
            ax1.set_yticklabels([])
            ax1.set_ylim([
                    -1*np.max(np.abs(sta_[cmp]["tr_results"][n_filt])),
                    np.max(np.abs(sta_[cmp]["tr_results"][n_filt]))
            ])
            for n_kurt, kurt in enumerate(KURT_WINS):
                ax3.plot(sta_[cmp]["times"],
                         sta_[cmp]["f3_results"][n_filt][n_kurt],
                         linewidth=0.7, linestyle="--")
        p = sta["best_obs"][0]
        for n in range(0, len(sta["best_obs"])):
            p = sta["best_obs"][n]
            if n == 0:
                phase = "P"
            else:
                phase = "S"
            if p:
                if p[2] == tr.stats.channel:
                    ax1.text(p[0]-orig_time, ax1.get_ylim()[1],
                             "{:1g}".format(p[3]),
                             horizontalalignment="center", bbox=dict(
                                     facecolor='none', edgecolor='red'))
                if cmp in cmps[phase]:
                    ax1.axvline(p[0]-orig_time, color=co[n], linestyle="--",
                                zorder=10,
                                label="Best Kurtosis {:}-pick".format(phase))
                    ax2.axvline(p[0]-orig_time, color=co[n], linestyle="--",
                                zorder=10)
                    ax3.axvline(p[0]-orig_time, color=co[n], linestyle="--",
                                zorder=10)
                    ax4.axvline(p[0]-orig_time, color=co[n], linestyle="--",
                                zorder=10)
                    ax1.set_xlim([
                            p[0]-orig_time-PLOT_WIND_PICK_DIFF,
                            p[0]-orig_time+PLOT_WIND_PICK_DIFF
                                  ])
                    ax2.set_xlim([
                            p[0]-orig_time-PLOT_WIND_PICK_DIFF,
                            p[0]-orig_time+PLOT_WIND_PICK_DIFF
                                  ])
                    ax3.set_xlim([
                            p[0]-orig_time-PLOT_WIND_PICK_DIFF,
                            p[0]-orig_time+PLOT_WIND_PICK_DIFF
                                  ])
                    ax4.set_xlim([
                            p[0]-orig_time-PLOT_WIND_PICK_DIFF,
                            p[0]-orig_time+PLOT_WIND_PICK_DIFF
                                  ])

        for phase in ["P", "S"]:
            if cmp in cmps[phase]:
                ax1.axvline(sta["pred"][phase]-orig_time,
                            label="Theor. {:}-arrival".format(phase),
                            color=ct[nphs], linestyle="--")
        for n_sw, sw in enumerate(CF_MEAN_SMOOTH_WIND):
            ax4.plot(sta_[cmp]["times"], sta_[cmp]["f4_all"][n_sw],
                     label="F$_4$ {:}".format(sw))
        ax4_ = ax3.twinx()
        ax3.plot(sta_[cmp]["times"], sta_[cmp]["f3_mean_smooth"][0],
                 color="red", linewidth=2, label="Smoothed mean CF")
        ax1.text(0.04, 0.96, "Waveforms {:}".format(tr.stats.channel),
                 horizontalalignment="left", verticalalignment="top",
                 transform=ax1.transAxes,
                 bbox=dict(facecolor="white", edgecolor="black",
                           boxstyle="round", alpha=0.5))
        ax2.text(0.04, 0.96, "Raw mean Kurtosis {:}".format(tr.stats.channel),
                 horizontalalignment="left", verticalalignment="top",
                 transform=ax2.transAxes, bbox=dict(
                        facecolor="white", edgecolor="black",
                        boxstyle="round", alpha=0.5))
        ax3.text(0.04, 0.96, "Kurtosis CF {:}".format(tr.stats.channel),
                 horizontalalignment="left", verticalalignment="top",
                 transform=ax3.transAxes, bbox=dict(
                        facecolor="white", edgecolor="black",
                        boxstyle="round", alpha=0.5))
        ax4.text(0.04, 0.96, "F$_4$ transform {:}".format(tr.stats.channel),
                 horizontalalignment="left", verticalalignment="top",
                 transform=ax4.transAxes, bbox=dict(
                        facecolor="white", edgecolor="black",
                        boxstyle="round", alpha=0.5))
        plt.suptitle("{:}.{:}".format(
                st_sta[n_tr].stats.network, st_sta[n_tr].stats.station),
                fontsize=16)
        ax1.grid()
        ax1.legend(loc="lower left")
        ax2.grid()
        ax2.legend(loc="lower left")
        ax3.grid()
        ax3.legend()
        ax4.grid()
        h1, l1 = ax4.get_legend_handles_labels()
        h2, l2 = ax4_.get_legend_handles_labels()
        ax4.legend(h1+h2, l1+l2, loc="upper right")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.subplots_adjust(wspace=0.1)
    pdf.savefig(fig)
    # plt.show()
    plt.cla()
    plt.clf()


def write_evt(evt, ev_id):
    """
    Write refined event with new picks to stationxml file
    """
    print("Writing xml file")
    print("../../refined_events/{:}.xml".format(ev_id))
    evt.write("../../refined_events/{:}.xml".format(ev_id), format="QUAKEML")
    #time = orig.time.isoformat()
    #orig = evt.preferred_origin()
    #lon = orig.longitude
    #lat = orig.latitude
    #depth = orig.depth/1000
    #mag = evt.preferred_magnitude().mag
    #rms = orig.quality.standard_error
    #nphs = orig.quality.associated_phase_count
    #gap = orig.quality.azimuthal_gap
    #w = open("refined_events.dat", "a")
    #w.write("{:} {:7.4f} {:7.4f} {:5.2f} {:5.2f} {:5.2f} {:3g} {:3g}\n".format(
    #    time, lat, lon, depth, mag, rms, nphs, gap))
    #w.close()

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


# Main function
# Read in catalogues
if __name__ == '__main__':
    initial_cat = read_events(lassie_cat_file)
    
    # Read in station locations from file
    sta_list = [[l.split()[1], float(l.split()[3]), float(l.split()[4])] for l in
                open("NLLOC_run/run.in", "r") if l.split()[0] == "GTSRCE"]
    sta_locs = {sta[0]: {"lat": sta[1], "lon": sta[2]} for sta in sta_list}
    if DO_PLOT_1 is True or DO_PLOT_2 is True:
        fig = plt.figure(figsize=(18, 10))
    else:
        fig = []
    
    if FORCE_RECALC is True:
        filelist = glob.glob(os.path.join("refined_events", "*.xml"))
        for f in filelist:
            os.remove(f)
    
    #nproc = multiprocessing.cpu_count()
    nproc = 4
    # Get events for which data currently is available
    cat_filter = Catalog()
    for n, event in enumerate(initial_cat):
        e_id = event.event_descriptions[0].text
        if e_id in os.listdir("{:}/".format(DIR_TO_EVENTDIRS)):  # and n == 3):
            if (FORCE_RECALC is False and
                    os.path.exists("refined_events/{:}.xml".format(e_id))):
                    print("Already have this evening ... skipping")
            else:
                cat_filter.append(event)
    # Split catalogue across multiple processes and process in parallel
    cat_split = [i for i in mit.divide(nproc, cat_filter)]
    #process_events(cat_split[7], 7)
    processes = []
    pool = multiprocessing.Pool(processes=nproc)
    print("hello")
    a = pool.starmap(process_events, product(cat_split, range(nproc)))

    logger.debug(a)
    #for i in range(nproc):
#    print("Running process", i)
#    mp = multiprocessing.Process(target=process_events,
#                                 args=(cat_split[i], i))
#    mp.start()
#    processes.append(mp)
#
#for p in processes:
#    p.join()
