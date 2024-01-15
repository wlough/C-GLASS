import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import proj3d

n_datapoints = 1000
n_vrts = 642
n_adj_max = 10

data_pos = np.fromfile(
    "/home/shane/projects/C-GLASS/test_membrane_vrt_positions.file",
    dtype=np.double,
).reshape(n_datapoints, n_vrts, 3)
data_adj = np.fromfile(
    "/home/shane/projects/C-GLASS/test_membrane_vrt_adjacency.file",
    dtype=np.int32,
).reshape(n_datapoints, n_vrts, n_adj_max)

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection="3d")

for i in range(0, n_datapoints, 10):
    ax.clear()
    x = data_pos[i, :, 0]
    y = data_pos[i, :, 1]
    z = data_pos[i, :, 2]
    ax.scatter(x, y, z)
    # TODO vectorize this
    for j in range(0, n_vrts):
        x0 = data_pos[i, j, 0]
        y0 = data_pos[i, j, 1]
        z0 = data_pos[i, j, 2]
        for k in range(0, n_adj_max):
            j_neighb = data_adj[i, j, k]
            if j_neighb == -1:
                continue
            x1 = data_pos[i, j_neighb, 0]
            y1 = data_pos[i, j_neighb, 1]
            z1 = data_pos[i, j_neighb, 2]
            ax.plot([x0, x1], [y0, y1], [z0, z1])
            ax.set_xlim(-40, 40)
            ax.set_ylim(-40, 40)
            ax.set_zlim(-40, 40)
    plt.pause(0.05)
plt.show()
