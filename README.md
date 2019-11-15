Help on module AutoQuakePycker:

NAME
    AutoQuakePycker

DESCRIPTION
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

FUNCTIONS
    CF_kurtosis(kurt_win_size, tr)
        This is the Kurtosis function of Baillard et al. (2014, SRL)
        (Eq. 10). Sets values outside of window to NaN.
    
    add_picks(sta, e, orig_time, st)
        Add picks made to ObsPy event structure.
        Args:
            sta (dict): variable station & picking data.
            e (ObsPy event object): event.
            orig_time (ObsPy UTCDateTime object): event origin time.
            st (Obspy stream object): stream of ObsPy traces containing data.
        Returns:
            e (ObsPy event object): event containg added picks.
    
    compute_magnitude(evt, st)
        Compute local magnitude for event. Uses the same approach as Bie et al.
        (2019, SRL), although we don't bother with checking signal-to-noise ratio
        this time.
        Args:
    
    compute_noise_levels(tr)
        Get simple estimate of trace noise level using STA/LTA function - similar
        to thhat described by Baillard et al. (2014, BSSA).
    
    copytree(src, dst, symlinks=False, ignore=None)
    
    do_smooth(d, WT, sample_rate)
        Simple windowing smoothing function.
        (takes data to right of smoothing point).
    
    get_best_picks(sta, max_wind, orig_time, comps, sptype)
        Define best Picks
        Choose whether to use fixed min S-P time or variable S-P time depending on
        P-arrival time and Vp/Vs ratio (Wadati method)
        This is quite messy. Could benefit from some improvements.
    
    get_theoretical_tt(event)
        Find theoretical arrival times from NLLoc Time2EQ (containing 1 event).
    
    kurt_transform_f2(f1, kurt_win_size, tr)
        Apply transformations as per Baillard (BSSA).
        Remove neg slope of the CF from the sliding kurtosis fonction.
    
    kurt_transform_f3(f2, kurt_win_size, tr)
        Remove the linear trend of F2.
        Eq. 12
    
    kurt_transform_f4(f3, kurt_win_size, tr)
        Find greatest minima correspond to the greatest onset strengths.
    
    plot_record_section(st_sta)
        TODO
    
    plot_sta_results(sta, st_sta, runtype, orig_time, pdf, fig, cmps)
        Plot picking results for each station
    
    process_events(cat_data, n_run)
        This is the main function that runs a given set of events and does the
        multi-stage refinement.
    
    process_station(sta, st_sta, dist, ev_dict, evt, orig_time, phase, cmps)
        Get ready to process data. Sets up emtpy arrays.
        Checks noise level of data.
    
    process_trace(n_tr, tr, sta, orig_time, cmps)
        Process picking on each trace parsed by process_station.
    
    refine_events(st, stations_dist, cmps, max_pick_diff, evt, do_plot, n_refine, fig, sptype, path_to_figs, ev_dict, ev_id, rot=False)
        This is one of the main functions of this code. It processes a refinement
        attempt for each event.
    
    relocate_nlloc(evt_refined)
        Carry out relocation using newly picked events in NonLinLoc.
    
    rotate(st, evt)
        Rotate raw 3-component data to ZRT coordinate system.
        Throws up flag if rotation cannot be done.
    
    write_evt(evt, ev_id)
        Write refined event with new picks to stationxml file

DATA
    CF_MEAN_SMOOTH_WIND = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0....
    CMPS_REFINE_1 = {'P': ['Z', 'H'], 'S': ['1', '2', 'N', 'E']}
    CMPS_REFINE_2 = {'P': ['Z', 'H'], 'S': ['T', 'R']}
    DIR_TO_EVENTDIRS = '../events/events'
    DO_PLOT_1 = True
    DO_PLOT_2 = True
    EvaluationMode = Enum(["manual", "automatic"])
    FILT_WINS = {'P': [[4, 12], [10, 18]], 'S': [[2, 6], [3, 8]]}
    FORCE_RECALC = False
    INTEGRATE_S = True
    KURT2WGHT = {'P': [6, 4.8, 4, 2.9, 0], 'S': [6, 4.8, 4, 2.9, 0]}
    KURT_WINS = [2, 4, 6]
    MAG_HIGHPASS = 2
    MAX_DIST = 275
    MAX_DIST_ML = 400.0
    MAX_DRIFT = 0.12
    MAX_NOISE_LEVEL = 50.0
    MAX_PICK_DIFF_REFINE1 = 12.0
    MAX_PICK_DIFF_REFINE2 = 3.5
    MAX_WINDOW = 150.0
    MIN_SP_TIME = 2.5
    ML_CAL = [1.0, 0.00301, 3.0]
    PLOT_WIND_PICK_DIFF = 4.0
    SNR_FREQ = [4, 8]
    SNR_SMOOTH_WIN = 1.0
    SNR_WIN = [5, 1]
    STA_BLACKLIST = ['BAUV', 'TRNT']
    T_ERRORS = [0.1, 0.2, 0.6, 1.2, 99.9]
    V_MIN = 3.0
    lassie_cat_file = '../new_event_detections.xml'
    logger = <Logger multiprocessing (DEBUG)>

FILE
    /Users/sph1r17/RESEARCH_PROJECTS/VOILA/2019_detection_picks/AutoQuakePycker/AutoQuakePycker.py


