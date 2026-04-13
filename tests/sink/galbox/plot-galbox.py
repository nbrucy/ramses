import osyris
import matplotlib.pyplot as plt
import visu_ramses


data = osyris.RamsesDataset(2, ".").load()

center = osyris.Vector(x=0.5, y=0.5, z=0.5, unit="kpc")

osyris.map(
    data["mesh"].layer("density"),
    data["sink"].layer("position", mode="scatter"),
    norm="log",
    dx=1000 * osyris.units("pc"),
    #dz=1000 * osyris.units("pc"),
    direction="z",
    origin=center,
    cmap="inferno",
)

plt.savefig('galbox.pdf',bbox_inches='tight')

dat_to_cmp = {}
for key in data["mesh"]:
    dat_to_cmp[key] = data["mesh"][key].values.flatten()
for key in data["sink"]:
    dat_to_cmp[key] = data["sink"][key].values.flatten()


# Check results against reference solution
visu_ramses.check_solution(dat_to_cmp,'galbox', threshold=1e-30, overwrite=False)
