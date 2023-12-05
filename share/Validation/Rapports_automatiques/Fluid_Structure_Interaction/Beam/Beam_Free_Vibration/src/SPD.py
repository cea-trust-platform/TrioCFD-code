import numpy as np


# Signal
t= np.loadtxt('Poutre_Displacement1D.out', unpack=True, usecols=[0])
x = np.loadtxt('Poutre_Displacement1D.out', unpack=True, usecols=[3])
max= np.max(x);
signal = x/max  # signal array

# Compute the FFT
fft = np.fft.fft(signal)

# Compute the power spectrum
power_spectrum = np.abs(fft)**2

# Find the index of the maximum value in the power spectrum
max_index = np.argmax(power_spectrum)

# Compute the frequency corresponding to the maximum value in the power spectrum
sampling_freq = 1 / (t[1] - t[0])  # sampling frequency
freq = max_index * sampling_freq / len(signal)

freq= round(freq,2)
fichier = open("Freq.txt", "w")       
fichier.write(str(freq))        
fichier.close()


