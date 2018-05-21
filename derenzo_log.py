import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpp

# Create figure with aspect ration = 1.0
fig = plt.figure(figsize=(6, 6))
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

# Create the outer cylinder of the phantom
R = 50.0
derenzo_cyl = mpp.Circle((0, 0), radius=R, color="gray", alpha=0.3)
ax.add_patch(derenzo_cyl)
ax.set_xlim((-1.2*R, 1.2*R))
ax.set_ylim((-1.2*R, 1.2*R))

# Partition phantom into 6 sections
hpts = np.array([(-R, 0), (R, 0)])
plt.plot(*hpts.T, color='k')
th = -60. * (np.pi/180)
rot = np.array([(np.cos(th), -1*np.sin(th)),
                (np.sin(th), np.cos(th))])
hpts = np.dot(hpts, rot)
plt.plot(*hpts.T, color='k')
hpts = np.dot(hpts, rot)
plt.plot(*hpts.T, color='k')

# For a given feature size/separation, compute the number of rows & total 
# number of wells in the section
section_offset = 0.10 * R   # Well region offset from the true section bdries
well_sep = 1.0
r = well_sep / 2.0          # Feature size related to separation, for now
h_section = R - (2 * section_offset + 2 * r)
row_height = well_sep * np.sqrt(3)
num_rows = int(np.floor(h_section / row_height))

# Compute well locations given the well separation and number of rows
num_wells = np.sum(1 + np.arange(num_rows))
xs, ys = [], []
for i in range(num_rows):
    rn = i + 1
    for x in np.arange(-rn, rn, 2) + 1:
        xs.append(x * well_sep)
        ys.append(-(section_offset + row_height * rn))
locs = np.vstack((xs, ys)).T

# Plot wells
for xy in locs:
    cyl = mpp.Circle(xy, radius=r, color="green", alpha=0.5)
    ax.add_patch(cyl)

plt.show()
