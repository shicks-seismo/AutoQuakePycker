#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Various wrapping functions for integrating AutoQuakePycker with NonLinLoc

@author: Stephen Hicks
"""
import os
import subprocess


def relocate(evt_refined):
    """
    Carry out relocation using newly picked events in NonLinLoc

    Parameters
    ----------
    evt_refined : ObsPy event object
    tr : ObsPy trace object

    Returns
    ----------
    ev_new : ObsPy event object containing relocated event
    """
    from obspy import read_events
    if len(evt_refined.picks) <= 5:
        ev_new = []
    else:
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


def get_theoretical_tt(event):
    """
    Find theoretical arrival times from NLLoc Time2EQ (containing 1 event).

    Parameters
    ----------
    event : ObsPy event object

    Returns
    ----------
    e : ObsPy event object containing theoretical arrival times
    """
    from obspy.core.event import Pick, WaveformStreamID

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
