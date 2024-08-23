import scipy.io
import matplotlib.pyplot as plt
import compSensPack as csp

if __name__ == "__main__":

    
    # I put this here because you have to decide kron fact based on
    # the signal length, but for testing reasons I prefer to do the 
    # opposite: choose signal length based on kron factor
    
    n_block = 16 # block size, deafult is 16
    kron_factor = n_block * 2  # kron factor default is 32


    ## DATA
    # ----------------------------------------------------------------
    #  upload data
    data = scipy.io.loadmat('100m.mat')
    # retrieve the key to a string
    key = list(data.keys())[0]
    # retrieve the values
    signal = data[key][0,:]  # [0 or 1, 0:650000] s.t. first dim: (0 is MLII, 1 is V5)
    num = 10  # how many kronecker blocks will be long the original signals
    temp = n_block * kron_factor  # this is equal to the length of a kronecker block
    start = int(temp * 0)  # choose where to start our signal in the record (signal is a piece of the record)
    end = int(start + temp * num)
    signal = signal[start:end]  # comment out to use the whole signal
    # print the shape of the signal
    print(f'Signal shape: {signal.shape}') 
    # plot the signal
    plt.figure(figsize=(10, 6))
    plt.plot(signal)
    plt.title('ECG Signal')
    plt.xlabel('Index')
    plt.ylabel('Amplitude')
    plt.grid(True)
    plt.show()    

    ## PARAMETERS
    # ----------------------------------------------------------------
    CR = 1/4  # compression ratio
    # NON-KRONECKER
    #n_block = 16 # block size`
    m_block = int(n_block * CR) # compressed block size
    # KRONECKER
    #kron_factor = 32  # kron factor
    n_block_kron = n_block * kron_factor  # kron block size
    m_block_kron = int(n_block_kron * CR) # compressed kron block size
    


    ### SAMPLING PHASE 
    # ----------------------------------------------------------------
    # MEASUREMENT MATRIX
    Phi = csp.generate_DBDD_matrix(m_block, n_block)
    #Phi = generate_random_matrix(m_block, n_block, matrix_type='gaussian')
    #Phi = generate_random_matrix(m_block, n_block, matrix_type='scaled_binary')
    #Phi = generate_random_matrix(m_block, n_block, matrix_type='unscaled_binary')

    # COMPRESS THE SIGNAL
    Y = csp.compressSignal(signal, Phi)



    ## RECOVERY PHASE
    # ----------------------------------------------------------------
    # SL0 PARAMETERS
    sigma_off = 0.001  # offset for sigma
    mu_0 = 2  # scaling factor for mu
    sigma_decrease_factor = 0.5  # factor for decreasing sigma
    L = 3  # number of iterations for the inner loop
    if sigma_off > 0:
        sigma_min = sigma_off * 4  # minimum value of sigma
    else:
        sigma_min = 0.00001  # minimum value of sigmZ
    # RECOVERY non-KRONECKER
    recovered_signal = csp.non_kron_recovery(Y=Y, sigma_min=sigma_min, Phi=Phi, sigma_decrease_factor=sigma_decrease_factor,
                                            mu_0=mu_0, L=L, showProgress=False)
    # RECOVERY KRONECKER
    recovered_signal_kron = csp.kron_recovery(Y=Y, sigma_min=sigma_min, Phi=Phi, sigma_decrease_factor=sigma_decrease_factor,
                                            mu_0=mu_0, L=L, showProgress=False, kron_factor=kron_factor)
    

    ## EVALUATION
    # ----------------------------------------------------------------
    # Plot and SNR
    csp.plot_signals(signal, recovered_signal, original_name="Original Signal", reconstructed_name="Recovered Signal")
    csp.plot_signals(signal, recovered_signal_kron, original_name="Original Signal", reconstructed_name="Recovered Signal (Kronecker)")
