import matplotlib.pyplot as plt
import numpy as np

no_of_elements = (3425, 17450, 41052, 89523, 136026, 177779)
tip_temp = (47.6741, 47.6407, 47.0462, 47.1739, 47.1735, 47.1199)
bulk_temp = (44.3402, 43.6827, 43.1146, 43.5092, 43.5981, 43.4977)

com_time = (22.7, 254.2, 274.5, 674.5, 1269.9, 2782.2)
x = (136026, 136026)

fig, ax1 = plt.subplots()

ax2 = ax1.twinx()
ax1.plot(no_of_elements, tip_temp, marker = '+', label="Tip point", color='C0')
ax1.plot(no_of_elements, bulk_temp, marker = 'o', label="Lower corner point", color='C0')
ax1.plot((136026, 136026),(42.5,49), ls = '--', color='k',  lw = 0.5)
ax2.plot(no_of_elements, com_time, ls = '--', label="Time", color='C1')
ax1.set_xlabel('Number of Elements', fontsize = 10)
ax1.set_ylabel('Temperatures ($^\circ$C)', color='C0', fontsize = 10)
ax2.set_ylabel('Computation Time (sec)', color='C1', fontsize = 10)
ax1.set_xlim(left = 0)
ax1.set_ylim(42.5,49)
ax2.set_ylim(0,3000)
ax1.tick_params(labelsize = 10)
ax2.tick_params(labelsize = 10)
ax1.legend(loc = "upper left", fontsize = 10)
ax2.legend(loc = "upper right", fontsize = 10)
plt.show()
