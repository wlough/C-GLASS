import numpy as np
import matplotlib.pyplot as plt

n_datapoints = 250

data = np.fromfile(
    "/home/shane/projects/C-GLASS/test_membrane_forces.file", dtype=np.double
).reshape(n_datapoints, 4)

t = np.arange(0, n_datapoints)

fig, ax = plt.subplots()
ax.plot(t, data[:, 0], label="teth")
ax.plot(t, data[:, 1], label="bend")
ax.plot(t, data[:, 2], label="area")
ax.plot(t, data[:, 3], label="vol")
ax.set_xlim([0, n_datapoints])
# ax.set_ylim([-5, 50])
plt.legend(loc="upper right")

ax.set(xlabel="time (s)", ylabel="force (AU)")
ax.grid()

fig.savefig("test_1x.png", dpi=300)
plt.show()
