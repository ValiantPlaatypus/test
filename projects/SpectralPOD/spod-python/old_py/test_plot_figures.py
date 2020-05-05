import matplotlib.pyplot as plt
import numpy as np
import multiprocessing

data = np.random.random((100, 1000))
x = list(range(100))

fig = plt.figure()
ax = fig.add_subplot(111)

def save_plot(i):
    ax.scatter(x, data[:,i])
    fig.savefig("%d.png" % i, dpi=100)
    ax.cla()

p = multiprocessing.Pool(4)
p.map(save_plot, range(data.shape[1]))
