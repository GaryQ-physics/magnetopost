from events import eventlist
from matplotlib import pyplot as plt

event_list = eventlist()

plt.figure()
plt.hist(event_list[:,7], bins=range(25), density=True, color='k')
plt.xticks(range(0, 27, 3))
plt.xlabel('MLT')
plt.xlim([0, 24])
plt.grid()
plt.show()

plt.figure()
plt.hist(event_list[:,5], bins=range(45, 90, 5), density=True, color='k')
plt.xticks(range(45, 95, 5))
plt.xlabel('Magnetic Latitude [degrees]')
plt.xlim([45, 90])
plt.grid()
plt.show()


import numpy as np
(u, idx, n) = np.unique(event_list[:, 0:5], axis=0, return_index=True, return_counts=True)
event_list = event_list[idx, :]
