import numpy as np

def round_up_to_power_of_two(n):
    p = 1
    while p < n:
        p += p
    return p


def correlation(ts1, ts2=None):
    """
    Correlation function
    :param ts1: time series (scalar or vector)
    :type ts1: np.ndarray
    :param ts2: second time series (scalar or vector), if None, assume ts2=ts1
    :type ts2: np.ndarray
    :returns: the correlation function of ts1 and ts2
    :rtype: np.ndarray
    """

    n = len(ts1)
    n2 = 2*round_up_to_power_of_two(n)

    fft_ts1 = np.fft.fft(ts1, n2, axis=0)

    if ts2 is None:
        # auto-correlation
        ts2 = ts1
        fft_ts2 = fft_ts1
    else:
        # cross-correlation
        assert ts2.shape == ts1.shape
        fft_ts2 = np.fft.fft(ts2, n2, axis=0)

    product = np.conjugate(fft_ts1)*fft_ts2
    assert len(product) == n2
    fft_corr = np.fft.ifft(product, n2, axis=0)

    if len(fft_corr.shape) == 1:
        # a scalar time series
        return fft_corr.real[:n] / np.arange(n, 0, -1)
    else:
        # a vector time series - implied dot product
        return np.add.reduce(fft_corr.real[:n], axis=1) / np.arange(n, 0, -1)
