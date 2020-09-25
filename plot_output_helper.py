import numpy as np
import matplotlib.pyplot as plt

md_data = np.loadtxt("w6_dynamics_production_full.md", skiprows=7)

print(np.mean(md_data[:,5]))

def running_average(data, window_width):
    cumsum_vec = np.cumsum(np.insert(data, 0, 0))
    return (cumsum_vec[window_width:] - cumsum_vec[:-window_width]) / window_width


plt.xlabel("Running Average over 10000 Steps")
plt.ylabel("Running Average Kinetic Temperature (K)")
plt.title(r'Cage Water Hexamer, dt=0.5fs')
running_average_temp = running_average(md_data[:,5], 1000)
plt.plot(range(len(running_average_temp)), running_average_temp)
plt.axhline(y=np.mean(md_data[:,5]), color='r', linestyle='-')
plt.show()
