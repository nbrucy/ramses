import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import osyris as visu_ramses


# Make figure
fig = plt.figure()
ratio = 0.75
sizex = 12.0
fig.set_size_inches(sizex,ratio*sizex)

ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)

# Load RAMSES output
data = visu_ramses.RamsesDataset(2)
data.load()

x = data["mesh"]["position"].x
y = data["mesh"]["position"].y
r = np.sqrt(x**2 + y**2)
rho = data["mesh"]["density"]
P = data["mesh"]["pressure"]

ux  = data["mesh"]["velocity"].x
uy  = data["mesh"]["velocity"].y

vr = ux*x + uy*y # radial velocity
vphi = -uy*x + ux*y # azimuthal velocity

# 2D maps

mesh = data["mesh"]

ind = np.argmax(mesh["density"])
center = mesh["position"][ind]

visu_ramses.map(mesh.layer("density"),
    norm="log",
    origin=center,
    direction="z",
    ax=ax1,
)

visu_ramses.map(mesh.layer("pressure"),
    norm="log",
    origin=center,
    direction="z",
    ax=ax2,
)

visu_ramses.map(
    mesh.layer("velocity"),
    mode="vec",
    origin=center,
    direction="z",
    color="k",
    ax=ax3,
)

visu_ramses.map(
    mesh.layer("grav_acceleration"),
    mode="vec",
    origin=center,
    direction="z",
    color="k",
    ax=ax4,
)




fig.subplots_adjust(wspace=0.35)
fig.savefig('alpha_disk_1body.png',bbox_inches='tight')

# Check results against reference solution
#visu_ramses.check_solution(data["data"], 'alpha_disk_1body', tolerance={"all":3.0e-06}, overwrite=True)
