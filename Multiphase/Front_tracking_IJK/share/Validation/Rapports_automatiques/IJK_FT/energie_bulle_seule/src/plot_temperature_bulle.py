import numpy as np
import matplotlib.pyplot as plt

temp = np.loadtxt('basic_temperature_bulles_Ti_0.out')
plt.figure()
plt.plot(temp[:,0], temp[:,1], label='Temperature')
plt.xlabel(r'Time (s)')
plt.ylabel(r'Temperature (K)')
plt.show()
