import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import osyris

fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(12, 8))

# Load the data
data = osyris.RamsesDataset(2).load()
data0 = osyris.RamsesDataset(1).load()

rd = 0.125 * osyris.units("cm") 


center = osyris.Vector(x=1, y=1, z=1, unit="cm")

mesh = data['mesh']
mesh0 = data0['mesh']

for m in [mesh0, mesh]:
    m["centered_position"] = m["position"] - center
    m["r"] = m["centered_position"].norm
    m["vr"] = m["velocity"].dot(m["centered_position"])/m["r"]

osyris.map(
    mesh0.layer('density', norm="log"),
    ax=ax[0, 0],
    cmap='viridis',
    direction='z',
    origin=center,
    dx=2*rd,
)


osyris.map(
    mesh0.layer('density', norm="log"),
    ax=ax[0, 1],
    cmap='viridis',
    direction='x',
    origin=center,
    dx=2*rd,
)


osyris.map(
    mesh0.layer('vr', norm="log"),
    ax=ax[0, 2],
    cmap='viridis',
    direction='z',
    origin=center,
    dx=2*rd,
)



osyris.map(
    mesh0.layer('velocity'),
    ax=ax[0, 1],
    cmap='viridis',
    direction='x',
    origin=center,
    dx=2*rd,
    norm="log",
)

osyris.map(
    mesh.layer('grav_acceleration'),
    ax=ax[1, 0],
    cmap='viridis',
    direction='z',
    origin=center,
    dx=2*rd,
    norm="log",
)



osyris.map(
    mesh.layer('grav_acceleration'),
    ax=ax[0, 1],
    cmap='viridis',
    direction='x',
    origin=center,
    dx=2*rd,
    norm="log",
)


osyris.map(
    mesh.layer('density', norm="log"),
    ax=ax[2, 0],
    cmap='viridis',
    direction='z',
    origin=center,
    dx=2*rd,
)




osyris.map(
    mesh.layer('density', norm="log"),
    ax=ax[2, 1],
    cmap='viridis',
    direction='x',
    origin=center,
    dx=2*rd,
)


osyris.map(
    mesh.layer('vr'),
    ax=ax[2, 2],
    cmap='viridis',
    direction='z',
    origin=center,
    dx=2*rd,
)

for axi in ax.flat:
    axi.set_aspect('equal')
    
    


fig.savefig('betadisk.png',bbox_inches='tight')
