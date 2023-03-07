#!/Volumes/A1/workdir/nicole/bin/python
import xarray as xr
import matplotlib.pyplot as plt
import os

dir='/home/nicole/workdir/SWA14/forcing_subset_ERA5/'
files=os.listdir(f'{dir}')
print(files)
for f in files:
    a=xr.open_dataset(f'{dir}{f}')
    print(f'{dir}{f}')
#    nz=len(a.zl.data)
    nx=len(a.longitude.data)
    ny=len(a.latitude.data)

    cmaps = {'t2m': 'seismic', 'msl':'PiYG', 'u10':'nipy_spectral', 'v10':'nipy_spectral',
            'huss': 'gist_earth', 'siconc': 'bwr', 'sst':'seismic', 'ssrd':'jet', 'strd':'jet',
            'trr':'bwr_r', None:'viridis'}

    for v in a.keys():
        print(v)
    
    try:
        a[v][0,0,:,:].plot(size=10,cmap=cmaps[v])
        plt.savefig(f"FORC_{f}.png")
        plt.clf()
    except:
        a[v][0,:,:].plot(size=10,cmap=cmaps[v])
        plt.savefig(f"FORC_{f}.png")
        plt.clf()

print('done')
quit()
