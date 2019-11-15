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
