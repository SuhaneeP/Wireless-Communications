import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm


def generate_BER_vs_EbN0_Th(EbN0dB_range):
    """Generate the theoretical BER array for a given signal and SNR range

    The values that fall in the given SNR range are converted to linear form and the
    formula of BER for (1*2) SIMO MRC is used to compute the BER value.

    Args:
        EbN0dB_range (range): range object containing the range of SNR
                              values used in simulations

    Returns:
        BER (list): list of BER values
    """
    BER = []
    for EbN0dB in EbN0dB_range:
        EbN0 = 10.0 ** (EbN0dB / 10.0)
        mu = np.sqrt(EbN0 / (EbN0 + 1))
        Pe = 0.25 * (2 - 3 * mu + mu ** 3)
        BER.append(Pe)

    return BER


def generate_noise(signal, EbN0):
    """Generate the noise array for a given signal and its SNR

    The array of randomly generated values is multiplied with the noise variance which
    is computed using the SNR value. The array thus obtained is the noise array.

    Args:
        signal (array): signal array
        EbN0 (float): SNR value

    Returns:
        noise (array): The noise array
    """
    noise = (1 / np.sqrt(2 * EbN0)) * (
        np.random.normal(size=(np.shape(signal)))
        + np.random.normal(size=(np.shape(signal))) * 1j
    )
    return noise


def generate_BER_vs_EbN0_Pr(EbN0dB_range):
    """Generate practical BER for a BPSK signal

    The randomly generated message is converted to its corresponding BPSK signal by
    changing all the 0s to -1. The resulting signal array denotes the phase shift
    property. The same signal is duplicated nRx times and the final signal array is
    generated. For wireless channel, we multiply the signal array with the h array
    which is an array of complex numbers belonging to the normal distribution. Noise
    of the same dimension is added to the signal. On the receiver end, we multiply the
    signal received on each antenna by their corresponding channel factor. The signals
    received are summed up and by checking the real parts of the sum array, the received
    message is obtained.

    Args:
        EbN0dB_range (range): range object containing the range of SNR
                              values used in simulations

    Returns:
        BER (list): list of BER arrays where each array contains BER values for a
                    (1*nRx) SIMO MRC system
    """
    msg = np.random.randint(low=0, high=2, size=int(1e6))
    signal = np.ones(len(msg))
    for i, j in enumerate(msg):
        if j == 0:
            signal[i] = -1
    x = []
    x.append(signal)
    signal = np.asarray(x)

    ber = []
    for nRx in range(1, 5):
        ber_nRx = []
        for EbN0dB in EbN0dB_range:
            EbN0 = 10 ** (EbN0dB / 10.0)
            signal_prime = np.repeat(signal, nRx, axis=0)
            n = generate_noise(signal_prime, EbN0)
            h = (1 / np.sqrt(2)) * (
                np.random.normal(size=(np.shape(signal_prime)))
                + 1j * np.random.normal(size=(np.shape(signal_prime)))
            )
            y = h * x + n
            d = np.sum(np.conjugate(h) * y, axis=0)
            received_msg = np.real(d) >= 0
            Pb_pr = np.count_nonzero(msg != received_msg) / len(msg)
            ber_nRx.append(Pb_pr)
        ber.append(ber_nRx)

    return ber


if __name__ == "__main__":
    EbN0dB_range = range(0, 36)
    ber_Th_nRx1 = [
        0.5 * (1 - np.sqrt(10.0 ** (EbN0dB / 10.0) / (10.0 ** (EbN0dB / 10.0) + 1)))
        for EbN0dB in EbN0dB_range
    ]
    ber_Th_nRx2 = generate_BER_vs_EbN0_Th(EbN0dB_range)
    ber_Pr = generate_BER_vs_EbN0_Pr(EbN0dB_range)
    plt.plot(EbN0dB_range, ber_Th_nRx1, "x-", label="nRx=1 Analytical")
    plt.plot(EbN0dB_range, ber_Pr[0], "o-", label="nRx=1 Simulation")
    plt.plot(EbN0dB_range, ber_Th_nRx2, "x-", label="nRx=2 Analytical")
    plt.plot(EbN0dB_range, ber_Pr[1], "o-", label="nRx=2 Simulation")
    # plt.plot(EbN0dB_range, ber_Pr[2], ".-", label="nRx=3 Simulation")
    plt.plot(EbN0dB_range, ber_Pr[3], "o-", label="nRx=4 Simulation")
    plt.xscale("linear")
    plt.xlabel("SNR (dB)")
    plt.ylabel("BER")
    plt.ylim(bottom=1e-5)
    plt.yscale("log")
    plt.legend()
    plt.grid()
    plt.show()
