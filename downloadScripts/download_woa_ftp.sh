#!/Volumes/A1/workdir/nicole/bin/python

import os

# example ftp link temperature: https://www.ncei.noaa.gov/thredds-ocean/fileServer/woa23/DATA/temperature/netcdf/decav91C0/0.25/woa23_decav91C0_t16_04.nc
# salinity https://www.ncei.noaa.gov/thredds-ocean/catalog/woa23/DATA/salinity/netcdf/decav91C0/0.25/catalog.html

ftp_link='https://www.ncei.noaa.gov/thredds-ocean/fileServer/'

str="<a href='catalog.html?dataset="
Fin=[] ; Fout=[]
with open('catalog.html.4', 'r') as f:
    for line in f.readlines():
        if str in line:
            fle=line.replace(str,'').replace('>','').replace('<','').replace('tt','').replace('/a','').replace('/td','').split("'")
            print(fle[0],fle[1].replace('/td',''))
            Fin.append(fle[0])
            Fout.append(fle[1].replace('/',''))

for fin, fout in zip(Fin, Fout):
        if not os.path.isfile(f'/home/nicole/workdir/WOA/original/{fout}'):
                print(f"{fout}") 
                DATA=f"{ftp_link}/{fin}"
                os.system(f"wget {DATA} -O /home/nicole/workdir/WOA/original/{fout}")
quit()
