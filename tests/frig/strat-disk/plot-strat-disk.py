# REMARK: it is normal the reference value for density is negative. This is the sum of the log of the cell densities.

import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import osyris
import visu_ramses
from scipy.interpolate import griddata


# Load ramses's data
dataset = osyris.RamsesDataset(2)
dataset.load()

mesh = dataset["mesh"]
sink = dataset["sink"]
kpc = osyris.units("kpc")
center = osyris.Vector(x=0.5, y=0.5, z=0.5, unit=kpc)
fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(12, 8), sharex=True, sharey=True)
for i, direction in enumerate(["x", "y", "z"]):
    osyris.map(
        mesh.layer("density", norm="log"),
        direction=direction,
        dx=1 *osyris.units("kpc"), 
        dz=1 *osyris.units("kpc"),
        origin=center, 
        cmap="inferno",
        ax=ax[0][i],
    )
    
for i, direction in enumerate(["x", "y", "z"]):
    osyris.map(
        mesh.layer("velocity", norm="linear"),
        direction=direction,
        dx=1 *osyris.units("kpc"), 
        dz=1 *osyris.units("kpc"),
        origin=center, 
        cmap="Greens",
        ax=ax[1][i],
    )
    
for i, direction in enumerate(["x", "y", "z"]):
    osyris.map(
        mesh.layer("pressure", norm="log"),
        direction=direction,
        dx=1 *osyris.units("kpc"), 
        dz=1 *osyris.units("kpc"),
        origin=center, 
        cmap="Reds",
        ax=ax[2][i],
    )
    
ax[0][0].set_aspect("equal")
    
fig.savefig('strat-disk.pdf',bbox_inches='tight')

# Check results against reference solution
data = visu_ramses.load_data(2)
visu_ramses.check_solution(dataset["data"],'strat-disk', threshold=1e-30, overwrite=True)
