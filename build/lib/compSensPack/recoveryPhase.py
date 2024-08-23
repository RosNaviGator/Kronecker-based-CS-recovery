import numpy as np

from .sparseDictionaries import generate_DCT_dictionary  # Importing from sparseDictionaries.py
from .reconstructionMethods import SL0  # Importing from reconstructionMethods.py


def non_kron_recovery(Y, sigma_min, Phi, sigma_decrease_factor=0.5, mu_0=2, L=3, showProgress=False):

    n_block = Phi.shape[1]  # length of ORIGINAL signal block
    Dict = generate_DCT_dictionary(n_block)  # Dictionary
    Theta = Phi @ Dict  # Theta
    Theta_pinv = np.linalg.pinv(Theta)  # More-Penrose pseudoinverse of Theta
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