#!/usr/bin/env python
"""
mainscript.py
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

from running_funcs import process_events

from obspy import read_events, Catalog
import os
import glob
import multiprocessing
import more_itertools as mit
import logging
from itertools import product
from munch import munchify
import yaml

# Main function
# Read in catalogues
if __name__ == '__main__':
    logger = multiprocessing.log_to_stderr(logging.DEBUG)
    cfg = munchify(yaml.safe_load("config.yaml"))
    initial_cat = read_events(cfg.input.lassie_cat_file)

    # Read in station locations from file
    sta_list = [[l.split()[1], float(l.split()[3]), float(l.split()[4])] for l in
                open("NLLOC_run/run.in", "r") if l.split()[0] == "GTSRCE"]
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
    #process_events(cat_split[7], 7)
    processes = []
    pool = multiprocessing.Pool(processes=nproc)
    print("hello")
    a = pool.starmap(process_events, product(cat_split, range(nproc), cfg,
                                             sta_locs))

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
