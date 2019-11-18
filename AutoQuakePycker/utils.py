#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Routines for plotting the results of AutoQuakePycker
"""
from obspy import read_inventory
import numpy as np


def rotate(st, evt):
    """
    Rotate raw 3-component data to ZRT coordinate system.
    Throws up flag if rotation cannot be done
    (e.g. due to not all components spanning the same time)

    Parameters
    ----------
    st : ObsPy stream containing waveform data
    evt : ObsPy event object containing epicentre information for rotation

    Returns
    ----------
    st_sta_rot : ObsPy event object containing rotated traces
    done : Boolean flag stating whether rotation was succesful or not
    """
    from obspy.geodetics.base import gps2dist_azimuth
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


def add_picks(sta, e, orig_time, st, cfg):
    """
    Add picks made to ObsPy event structure.

    Parameters
    ----------
        sta : variable station & picking data (dict)
        e : event (ObsPy event object)
        orig_time : event origin time (ObsPy UTCDateTime object)
        st : stream of ObsPy traces containing data  (Obspy stream object)

    Returns
    ----------
        e : event containg added picks (ObsPy event object).
    """
    from obspy.core.event import (Pick, WaveformStreamID, EvaluationMode,
                                  QuantityError)
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
                time_errors=QuantityError(
                        uncertainty=cfg.picking.T_ERRORS[p[3]]))
            e.picks.append(pick)
    return(e)


def compute_magnitude(evt, st, cfg):
    """
    Compute local magnitude for event. Uses the same approach as Bie et al.
    (2019, SRL), although we don't bother with checking signal-to-noise ratio
    this time.

    Parameters
    ----------
    evt : ObsPy event object
    st : Obspy event stream containing waveform data
    cfg : Attribute-style dictionary containing picker configuration

    Returns
    ----------
    evt : ObsPy event object containing computed magnitude
    """
    from scipy.stats import trim_mean
    from obspy.core.event import Amplitude, WaveformStreamID, Magnitude

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
        st_mag.filter("highpass", freq=cfg.mag_calc.MAG_HIGHPASS)
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
        if dist_km >= cfg.mag_cal.MAX_DIST_ML:
            continue
        Lg_end = new_origin.time + dist_km / cfg.mag_cac.lV_MIN + 30
        if Lg_end - new_origin.time >= cfg.mag_calc.MAX_WINDOW:
            Lg_end = new_origin.time + cfg.mag_calc.MAX_WINDOW
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
        ML = (np.log10(amp) + cfg.mag_calc.ML_CAL[0] * np.log10(hyp_dist/100) +
              cfg.mag_calc.ML_CAL[1]*(hyp_dist-100) + cfg.mag_calc.ML_CAL[2]
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


def compute_noise_levels(tr, cfg):
    """
    Get simple estimate of trace noise level using STA/LTA function - similar
    to thhat described by Baillard et al. (2014, BSSA).
    """
    from obspy.signal.trigger import classic_sta_lta
    tr_snr = tr.copy()
    tr_snr.filter("bandpass", freqmin=cfg.sig_noise.SNR_FREQ[0],
                  freqmax=cfg.sig_noise.SNR_FREQ[1])
    wa = int(cfg.sig_noise.SNR_WIN[1]*tr.stats.sampling_rate)
    wb = int(cfg.sig_noise.SNR_WIN[0]*tr.stats.sampling_rate)
    # Prevent failing due to en(data) < nlta error
    if len(tr_snr.data) < wa or len(tr_snr.data) < wb:
        noise_level = 100.0
        return noise_level
    snr = classic_sta_lta(tr_snr.data, wa, wb)
    snr_smooth = do_smooth(snr, cfg.sig_noise.SNR_SMOOTH_WIN,
                           tr.stats.sampling_rate)
    thresh_snr = np.nanmax(snr_smooth) * 0.4
    A = (snr_smooth - thresh_snr)
    A = A[np.where(A > 0)]
    if len(snr_smooth[wb:-wa]) == 0:  # In case zerodivision error
        noise_level = 9999.9
        return noise_level
    noise_level = (len(A) / len(snr_smooth[wb:-wa])) * 100
    return noise_level


def copytree(src, dst, symlinks=False, ignore=None):
    """
    Make NLLOC run directory for each process

    Parameters
    ----------
    src : Copy source location (str)
    dst : Copy desination location (dst)
    """
    import os
    import shutil
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            shutil.copytree(s, d, symlinks, ignore)
        else:
            shutil.copy2(s, d)


def do_smooth(d, WT, sample_rate):
    """
    Simple windowing smoothing function.
    (takes data to right of smoothing point).

    Parameters
    ----------
        d : Time series to smooth (NumPy array)
        WT : Smoothing window length in seconds (float)
        sample_rate: sampling rate of time series (float)

    Returns
    ----------
        d_smooth : smoothed time series (NumPy array).
    """
    d_smooth = np.zeros(len(d))
    Wt = int(np.ceil(sample_rate*WT))
    for i in range(len(d)-Wt):
        d_smooth[i] = np.mean(d[i: i+Wt])
    d_smooth[0:Wt+100] = np.nan  # +100 removes "edge effects" at start of f4
    return(d_smooth)


def write_evt(evt, ev_id):
    """
    Write refined event with new picks to stationxml file

    Parameters
    ----------
        evt : Obspy event object
        ev_id : Event ID (str)
    """
    print("Writing xml file")
    print("../../refined_events/{:}.xml".format(ev_id))
    evt.write("../../refined_events/{:}.xml".format(ev_id), format="QUAKEML")
