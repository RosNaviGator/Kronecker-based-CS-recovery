import numpy as np

def compressSignal(signal, Phi):
    """
    Performs the sampling phase of Compressed Sensing (CS) by compressing a signal using a provided measurement matrix.

    In the context of Compressed Sensing, this function takes an input signal and a measurement matrix (Phi) and 
    compresses the signal into smaller blocks. Each column of the output matrix `Y` represents the compressed 
    measurement of a block of the input signal. The measurement matrix `Phi` is used to sample the signal in a way 
    that allows for efficient reconstruction from fewer samples than traditional Nyquist sampling would require.

    Parameters
    ----------
    signal : numpy.ndarray
        The original signal to be compressed, represented as a 1-dimensional array.
    
    Phi : numpy.ndarray
        The measurement matrix used for compressed sensing, of shape (M, N), where:
        - M is the number of measurements (rows of the compressed signal).
        - N is the length of each block of the original signal to be sampled.

    Returns
    -------
    Y : numpy.ndarray
        The compressed signal matrix of shape (M, BLOCK_NUM), where BLOCK_NUM is the 
        number of non-overlapping blocks created by dividing the length of the input signal by `N`.

    Notes
    -----
    - The function processes the input signal in non-overlapping blocks of length `N`. If the length of the signal 
      is not an exact multiple of `N`, the remaining samples at the end of the signal that do not fit into a full 
      block are ignored.
    - The measurement matrix `Phi` should have a number of columns (N) less than or equal to the length of the signal.
      If `Phi` has more columns than the length of the signal block, a `ValueError` will be raised due to a shape mismatch.
    - This function is a key step in the sampling phase of Compressed Sensing, where a sparse signal is projected 
      into a lower-dimensional space using random or structured measurements.
    """
    
    # length of signal block
    N = Phi.shape[1]

    # length of compressed block
    M = Phi.shape[0]

    # number of blocks
    BLOCK_NUM = len(signal) // N

    # each column of Y is the compressed version of a block of signal
    Y = np.zeros((M, BLOCK_NUM))    
    for i in range(BLOCK_NUM):
        Y[:,i] = Phi @ signal[i*N:(i+1)*N]
    
    return Y
