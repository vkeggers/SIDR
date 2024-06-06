import matplotlib.pyplot as plt
import numpy as np

# Simulate data
data1 = np.random.normal(500, 30, 200)  # Data around 500
data2 = np.random.normal(20, 5, 200)    # Data around 20

fig, (ax, ax2) = plt.subplots(2, 1, sharex=True, figsize=(6,8))
fig.subplots_adjust(hspace=0.05)  # adjust space between axes

# Plot the same data on both axes
ax.boxplot([data1, data2], positions=[1, 2])
ax2.boxplot([data1, data2], positions=[1, 2])

# Zoom-in / limit the view to different portions of the data
ax.set_ylim(480, 520)  # outliers only
ax2.set_ylim(0, 40)    # most of the data

# Hide the spines between ax and ax2
ax.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax.xaxis.tick_top()
ax.tick_params(labeltop=False)  # don't put tick labels at the top
ax2.xaxis.tick_bottom()

# This function is used to create a diagonal lines in the plot to indicate the break
d = .015  # how big to make the diagonal lines in axes coordinates
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

# Display the plot
plt.show()
