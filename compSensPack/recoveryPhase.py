import numpy as np

from .sparseDictionaries import generate_DCT_dictionary  # Importing from sparseDictionaries.py
from .SL0 import SL0  # Importing from reconstructionMethods.py


def non_kron_recovery(Y, sigma_min, Phi, sigma_decrease_factor=0.5, mu_0=2, L=3, showProgress=False):
    """
    Performs the non-Kronecker recovery phase of compressed sensing using the SL0 (Smoothed L0 norm) algorithm.

    This function reconstructs the original signal from the compressed measurements `Y` using a dictionary
    generated from the Discrete Cosine Transform (DCT). The SL0 algorithm is employed to solve the sparse
    recovery problem by iteratively minimizing a smoothed approximation to the L0 norm.

    Parameters
    ----------
    Y : numpy.ndarray
        A 2D array of shape (M, BLOCK_NUM), where each column represents a compressed block of the original signal.
    
    sigma_min : float
        The minimum value for the smoothing parameter sigma in the SL0 algorithm. Controls the final precision
        of the sparse solution.
    
    Phi : numpy.ndarray
        The measurement matrix of shape (M, N) used for compressing the original signal blocks. `Phi` should have
        more columns than rows (M < N) for an underdetermined system.
    
    sigma_decrease_factor : float, optional (default=0.5)
        The factor by which sigma is multiplied in each iteration of the SL0 algorithm to decrease it gradually.
    
    mu_0 : float, optional (default=2)
        The initial step size parameter for the SL0 algorithm. Controls the speed of convergence.
    
    L : int, optional (default=3)
        The number of iterations for each fixed sigma value in the SL0 algorithm. Higher values may improve accuracy
        but increase computation time.
    
    showProgress : bool, optional (default=False)
        If True, displays progress information during the SL0 algorithm iterations.

    Returns
    -------
    x_hat : numpy.ndarray
        A 1D array containing the recovered signal of length equal to the original signal. The reconstructed
        signal is obtained by applying the inverse transformation using the DCT dictionary.

    Notes
    -----
    - This function requires an implementation of the SL0 (Smoothed L0 norm) algorithm to perform sparse signal recovery.
      Ensure that the `SL0` function is correctly implemented and available in the environment where this function is called.
    - The SL0 algorithm is known for its fast convergence in sparse signal recovery problems and is effective
      in many practical compressed sensing applications.
    - This function assumes that the input matrix `Y` has been properly compressed using the provided measurement
      matrix `Phi` and that the DCT basis is suitable for the signal's sparse representation.

    Author
    ------
    RosNaviGator (https://github.com/RosNaviGator), 2024
    """

    n_block = Phi.shape[1]  # length of ORIGINAL signal block
    Dict = generate_DCT_dictionary(n_block)  # Dictionary
    Theta = Phi @ Dict  # Theta
    Theta_pinv = np.linalg.pinv(Theta)  # More-Penrose pseudoinverse of Theta, needed for SL0
    BLOCK_NUM = Y.shape[1]  # number of blocks

    n_block = Theta.shape[1]  # length of ORIGINAL signal block

    # Initialize the recovered signal
    x_hat = np.zeros(BLOCK_NUM * n_block)  # Recovered signal will be as long as the original signal

    # Process each block
    for i in range(BLOCK_NUM):
        y = Y[:, i]
        s_block = SL0(y, Theta, sigma_min, sigma_decrease_factor, mu_0, L, Theta_pinv, showProgress)
        x_hat[i*n_block:(i+1)*n_block] = Dict @ s_block

    return x_hat


def kron_recovery(Y, sigma_min, Phi, kron_factor, sigma_decrease_factor=0.5, mu_0=2, L=3, showProgress=False):
    """
    Performs the Kronecker product-based recovery phase of compressed sensing using the SL0 (Smoothed L0 norm) algorithm.

    This function reconstructs the original signal from the compressed measurements `Y` using a Kronecker product
    of the measurement matrix `Phi` and an identity matrix. The Kronecker product structure allows for the recovery
    of signals that are sparse in a Kronecker dictionary basis. The SL0 algorithm is used to solve the sparse
    recovery problem iteratively.

    Parameters
    ----------
    Y : numpy.ndarray
        A 2D array where each block represents a compressed section of the original signal. The matrix `Y` should
        be structured in such a way that its columns align with the blocks defined by the Kronecker structure.
    
    sigma_min : float
        The minimum value for the smoothing parameter sigma in the SL0 algorithm. Controls the final precision
        of the sparse solution.
    
    Phi : numpy.ndarray
        The original measurement matrix of shape (M, N) used for compressing the signal blocks. This matrix
        is expanded using a Kronecker product for block-structured compressed sensing.
    
    kron_factor : int
        The factor used to determine the size of the Kronecker product. Determines the number of times the
        identity matrix is replicated in the Kronecker product to expand `Phi`.
    
    sigma_decrease_factor : float, optional (default=0.5)
        The factor by which sigma is multiplied in each iteration of the SL0 algorithm to decrease it gradually.
    
    mu_0 : float, optional (default=2)
        The initial step size parameter for the SL0 algorithm. Controls the speed of convergence.
    
    L : int, optional (default=3)
        The number of iterations for each fixed sigma value in the SL0 algorithm. Higher values may improve accuracy
        but increase computation time.
    
    showProgress : bool, optional (default=False)
        If True, displays progress information during the SL0 algorithm iterations.

    Returns
    -------
    x_hat_kron : numpy.ndarray
        A 1D array containing the recovered signal of length equal to the original signal. The reconstructed
        signal is obtained by applying the inverse transformation using the Kronecker-structured dictionary.

    Notes
    -----
    - This function requires an implementation of the SL0 (Smoothed L0 norm) algorithm to perform sparse signal recovery.
      Ensure that the `SL0` function is correctly implemented and available in the environment where this function is called.
    - The Kronecker product structure is useful in compressed sensing when dealing with signals that exhibit
      sparsity across multiple dimensions or when the measurement process has a block structure.
    - The SL0 algorithm is used due to its efficiency and effectiveness in solving L0 minimization problems
      in compressed sensing, where the goal is to recover the sparsest solution possible.
    - Ensure that the input `Y` and `Phi` are correctly aligned and compatible with the Kronecker product structure
      and the sparse recovery assumptions.

    Author
    ------
    RosNaviGator (https://github.com/RosNaviGator), 2024
    """
    

    Phi_kron = np.kron(np.eye(kron_factor), Phi)  # KRONECKER product of Phi (kroncker measurement matrix)

    n_block = Phi.shape[1]  # length of ORIGINAL signal block
    #m_block = Phi.shape[0]  # length of compressed signal block
    n_block_kron = Phi_kron.shape[1]  # length of ORIGINAL signal KRONECKER-block
    #m_block_kron = Phi_kron.shape[0]  # length of compressed signal KRONECKER-block

    BLOCK_NUM = Y.shape[1]  # number of blocks
    N = n_block * BLOCK_NUM  # length of ORIGINAL signal
    KRON_BLOCK_NUM = N // n_block_kron  # number of KRONECKER blocks

    Dict_kron = generate_DCT_dictionary(n_block_kron)  # Dictionary (kron)
    Theta_kron = Phi_kron @ Dict_kron  # Theta (kron)
    Theta_kron_pinv = np.linalg.pinv(Theta_kron)  # More-Penrose pseudoinverse of Theta (kron)

    # Initialize the recovered signal
    x_hat_kron = np.zeros(N)

    # Process each block
    # DO NOT FORGET MATLAB IS COL-MAJOR, PY IS ROW-MAJOR, that darn 'F'... took my a day to figure out
    for i in range(KRON_BLOCK_NUM):
        y_kron = Y[:, i*kron_factor:(i+1)*kron_factor].reshape(-1, order='F')
        s_kron_block = SL0(y_kron, Theta_kron, sigma_min, sigma_decrease_factor, mu_0, L, Theta_kron_pinv, showProgress)
        x_hat_kron[n_block_kron*i:(i+1)*n_block_kron] = Dict_kron @ s_kron_block

    return x_hat_kron