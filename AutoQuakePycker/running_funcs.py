#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Main functions for processing events on different processors
"""

from AutoQuakePycker.utils import (rotate, do_smooth, add_picks, copytree,
                                   compute_magnitude, write_evt,
                                   compute_noise_levels)
from AutoQuakePycker.pickerfuncs import (
        CF_kurtosis, kurt_transform_f2, kurt_transform_f3, kurt_transform_f4,
        get_best_picks)
from AutoQuakePycker.nlloc_wrapping import relocate, get_theoretical_tt
from AutoQuakePycker.plotting import plot_sta_results
import numpy as np


def AutoQuakePycker_run():
    """
    Starts up AutoQuakePycker and sets-up processes
    """
    import glob
    from itertools import product
    import multiprocessing
    import logging
    import more_itertools as mit
    import munchify
    import os
    import yaml
    from obspy import read_events, Catalog

    logger = multiprocessing.log_to_stderr(logging.DEBUG)
    # Read in config file
    with open("config.yaml", "r") as ymlfile:
        cfg = munchify(yaml.safe_load(ymlfile))
    initial_cat = read_events(cfg.input.lassie_cat_file)

    # Read in station locations from file
    sta_list = [[l.split()[1], float(l.split()[3]), float(l.split()[4])] for l
                in open("NLLOC_run/run.in", "r") if l.split()[0] == "GTSRCE"]
    sta_locs = {sta[0]: {"lat": sta[1], "lon": sta[2]} for sta in sta_list}
    if cfg.output.FORCE_RECALC is True:
        filelist = glob.glob(os.path.join("refined_events", "*.xml"))
        for f in filelist:
            os.remove(f)
    if cfg.run.nproc == "auto":
        nproc = multiprocessing.cpu_count()
    else:
        nproc = cfg.run.nproc
    # Get events for which data currently is available
    cat_filter = Catalog()
    for n, event in enumerate(initial_cat):
        e_id = event.event_descriptions[0].text
        if e_id in os.listdir("{:}/".format(cfg.input.DIR_TO_EVENTDIRS)):
            if (cfg.output.FORCE_RECALC is False and
                    os.path.exists("refined_events/{:}.xml".format(e_id))):
                print("Already have this evening ... skipping")
            else:
                cat_filter.append(event)
    # Split catalogue across multiple processes and process in parallel
    cat_split = [i for i in mit.divide(nproc, cat_filter)]
    # process_events(cat_split[7], 7)
    pool = multiprocessing.Pool(processes=nproc)
    print("hello")
    a = pool.starmap(process_events, product(cat_split, range(nproc), cfg,
                                             sta_locs))

    logger.debug(a)


def process_events(cat_data, n_run, cfg, sta_locs):
    """
    This is the main function that runs a given set of events and does the
    multi-stage refinement.
    """
    import time
    import os
    import shutil
    import sys
    import logging
    from obspy import read
    from obspy.geodetics.base import gps2dist_azimuth
    import matplotlib.pyplot as plt

    if cfg.output.FORCE_RECALC is True:
        w = open("refined_events.dat", "w")
        w.close()
    if cfg.plotting.DO_PLOT_1 is True or cfg.plotting.DO_PLOT_2 is True:
        fig = plt.figure(figsize=(18, 10))
    else:
        fig = []
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
        st = read("../../{:}/{:}/MSEED/*.msd".format(
                cfg.input.DIR_TO_EVENTDIRS, ev_id), format="MSEED")
        print(n_run, ev_id)
        for n_tr, tr in enumerate(st):
            if st[n_tr].stats.sampling_rate > 40.0:
                try:
                    st[n_tr].resample(40)
                except ZeroDivisionError:
                    continue
        st1, st2, st_mag = [st.copy(), st.copy(), st.copy()]
        # Append distance to trace
        stations_data = sorted(set([tr.stats.station for tr in st
                                    if tr.stats.station not in
                                    cfg.sta_select.STA_BLACKLIST]))
        stations_dist = {sta_code: gps2dist_azimuth(
            sta_locs[sta_code]["lat"], sta_locs[sta_code]["lon"],
            orig_lat, orig_lon)[0] for sta_code in stations_data
            if gps2dist_azimuth(
                sta_locs[sta_code]["lat"], sta_locs[sta_code]["lon"],
                orig_lat, orig_lon)[0]/1000 <= cfg.sta_select.MAX_DIST}
        path_to_figs = "../../{:}/{:}/figs".format(
                cfg.input.DIR_TO_EVENTDIRS, ev_id)
        if not os.path.exists(path_to_figs):
            os.mkdir(path_to_figs)
        print("Doing first refinement")
        sys.stdout.flush()
        if ("R" in cfg.picking.CMPS_REFINE_1["S"] or
                "T" in cfg.picking.CMPS_REFINE_1["S"]):
            rot = True
        else:
            rot = False
        evt_refine_1, rms, found = refine_events(
                st1, stations_dist, cfg.picking.CMPS_REFINE_1,
                cfg.picking.MAX_PICK_DIFF_REFINE1, ev,
                cfg.ploting.DO_PLOT_1, 1, fig, "const", path_to_figs, ev_dict,
                ev_id, cfg, rot
                )
        if found is False:
            continue
        print("RMS = ", rms)
        sys.stdout.flush()
        prev_rms = rms
        print("Doing second refinement")
        sys.stdout.flush()
        if ("R" in cfg.picking.CMPS_REFINE_2["S"] or
                "T" in cfg.picking.CMPS_REFINE_2["S"]):
            rot = True
        else:
            rot = False
        evt_refine_2, rms, found = refine_events(
                st2, stations_dist, cfg.picking.CMPS_REFINE_2,
                cfg.picking.MAX_PICK_DIFF_REFINE2, evt_refine_1,
                cfg.plotting.DO_PLOT_2, 2, fig, "dist", path_to_figs, ev_dict,
                ev_id, rot
                )
        if found is False:
            continue
        print("RMS = ", rms)
        if rms > prev_rms * 1.25:
            print("RMS is significantly increasing (*25%) - skipping event")
            continue
        prev_rms = rms
        evt_refine_2 = compute_magnitude(evt_refine_2, st_mag, cfg)
        write_evt(evt_refine_2, ev_id)
        end = time.time()
        print("Time taken for event: {:3.1f} mins".format((end-start)/60))


def process_station(sta, st_sta, dist, ev_dict, evt, orig_time, phase, cmps,
                    cfg):
    """
    Get ready to process data. Sets up emtpy arrays.
    Checks noise level of data.
    """
    sta["ncha"] = len(st_sta)
    sta["lenD"] = int(st_sta[0].stats.npts)
    sta["samplerate"] = st_sta[0].stats.sampling_rate
    sta["nwinsamp"] = int(np.ceil(cfg.picking.KURT_WINS[0]*sta["samplerate"]))
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
        if sta["noise_levels"][n_tr] < cfg.sig_noise.MAX_NOISE_LEVEL:
            sta = process_trace(n_tr, tr, sta, orig_time, cmps, cfg)
    return(sta)


def process_trace(n_tr, tr, sta, orig_time, cmps, cfg):
    """
    Process picking on each trace parsed by process_station.
    """
    cmp = tr.stats.channel[2:3]
    sta[cmp] = {}
    sta[cmp]["times"] = tr.times(reftime=orig_time)

    sta[cmp]["tr_results"] = np.zeros(
            (len(cfg.picking.FILT_WINS["P"]), sta["lenD"])
            )
    sta[cmp]["f1_results"] = np.zeros(
            (len(cfg.picking.FILT_WINS["P"]), len(cfg.picking.KURT_WINS),
             sta["lenD"])
            )
    sta[cmp]["f1_mean"] = np.zeros(sta["lenD"])
    sta[cmp]["f3_results"] = np.zeros(
            (len(cfg.picking.FILT_WINS["P"]),
             len(cfg.picking.KURT_WINS), sta["lenD"])
            )
    sta[cmp]["f3_mean_smooth"] = np.zeros(
            (len(cfg.picking.CF_MEAN_SMOOTH_WIND), sta["lenD"])
            )
    sta[cmp]["f4_all"] = np.zeros((len(cfg.picking.CF_MEAN_SMOOTH_WIND),
                                   sta["lenD"]))
    sta[cmp]["f1_mean_smooth"] = np.zeros(sta["lenD"])
    # Get suitable filters (exclude those fully outside Nyquist freq.)
    for phase in ["P", "S"]:
        if cmp in cmps[phase]:
            sta["picks"]["poss_obs"][phase][cmp] = {}
            sta[cmp]["filtwins_check"] = [
                    filt_win for filt_win in cfg.picking.FILT_WINS[phase]
                    if filt_win[0] < sta["samplerate"] / 2
                    ]
        if cfg.picking.INTEGRATE_S is True:
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
        for n_kurt, kurt_win_s in enumerate(cfg.picking.KURT_WINS):
            f1 = CF_kurtosis(kurt_win_s, tr)
            sta[cmp]["f1_results"][n_filt, n_kurt] = f1  # Needed for weights
            f2 = kurt_transform_f2(f1, kurt_win_s, tr)
            f3 = kurt_transform_f3(f2, kurt_win_s, tr)

            sta[cmp]["f3_results"][n_filt, n_kurt] = f3
    sta[cmp]["f1_mean"] = np.nanmean(sta[cmp]["f1_results"], axis=0)[0]
    sta[cmp]["f1_mean_smooth"] = do_smooth(
            sta[cmp]["f1_mean"], cfg.picking.CF_MEAN_SMOOTH_WIND[0],
            tr.stats.sampling_rate
            )
    # ^ Throws up a warning first time due to NaN slices
    # Compute mean CF and final kurtosis transform
    f3_mean = np.nanmean(sta[cmp]["f3_results"], axis=0)[0]

    for nsm, smooth_wind in enumerate(cfg.picking.CF_MEAN_SMOOTH_WIND):
        sta[cmp]["f3_mean_smooth"][nsm] = do_smooth(
                f3_mean, smooth_wind, tr.stats.sampling_rate
                )
        f4 = kurt_transform_f4(sta[cmp]["f3_mean_smooth"][nsm],
                               np.max(cfg.picking.KURT_WINS), tr)
        sta[cmp]["f4_all"][nsm] = f4

        # Now pick (avoiding end and beginning of signal)
        # Pick the P-waves
        if cmp in cmps["P"]:
            sta["picks"]["poss_obs"]["P"][cmp][nsm] = []
            # Find points where Kurt<0 & doesn't look like S-wave
            p_cands = np.argwhere((f4 < 0.0))
            for idx in p_cands.tolist():
                kurt_wgt = np.min(np.where(np.array(
                        cfg.picking.KURT2WGHT["P"]
                        <= sta[cmp]["f1_mean_smooth"][idx])))
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
                kurt_wgt = np.min(np.where(np.array(cfg.picking.KURT2WGHT["S"]
                                           <= sta[cmp]["f1_mean_smooth"][idx]))
                                  )
                sta["picks"]["poss_obs"]["S"][cmp][nsm].append([
                        orig_time+sta[cmp]["times"][idx][0], f4[idx][0],
                        tr.stats.channel, kurt_wgt, idx,
                        sta[cmp]["times"][idx][0]
                        ])
    return(sta)


def refine_events(st, stations_dist, cmps, max_pick_diff, evt, do_plot,
                  n_refine, fig, sptype, path_to_figs, ev_dict, ev_id, cfg,
                  rot=False):
    """
    This is one of the main functions of this code. It processes a refinement
    attempt for each event.
    """
    import matplotlib
    from obspy import Stream
    from obspy.core.event import EventDescription

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
                        (sta["pred"][phase] - max_pick_diff -
                         np.max(cfg.picking.KURT_WINS) - 1),
                        (sta["pred"][phase] + max_pick_diff +
                         np.max(cfg.picking.KURT_WINS)))
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
                                     cmps, sptype, cfg)
            except:
                continue
            if sta["best_obs"][0]:
                evt = add_picks(sta, evt, orig_time, st_, cfg)
                if do_plot is True:
                    plot_sta_results(sta, st_sta, "refine{:}".format(
                            n_refine), orig_time, pdf, fig, cmps, cfg)

    if do_plot is True:
        pdf.close()
    evt_refine_out = relocate(evt)
    if evt_refine_out == []:
        print("Relocation failed - skipping event")
        found = False
        rms = 100.0
    else:
        evt_refine_out.event_descriptions = [(EventDescription(text=ev_id))]
        rms = evt_refine_out.preferred_origin().quality.standard_error
        found = True
    return(evt_refine_out, rms, found)
