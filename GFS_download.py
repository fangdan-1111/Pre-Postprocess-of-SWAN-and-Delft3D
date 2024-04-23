import xarray as xr
import datetime as dt
import numpy as np
import wget
import os

# setting
lat = (12, 42)
lon = (99, 142)
preci = 1
time = (0, 168, 6)
runs_start = [0]
Dir=r"\\wsl.localhost\Ubuntu-20.04\home\cxl\WRF_WorkStation\WPS-4.5\runcase"

# initialize
runs = runs_start+np.arange(*time)
range1 = lambda start, end, step=1: np.arange(start, end + step, step)
date1 = dt.datetime.now() - dt.timedelta(days=1)
date1 = date1.replace(hour=0,minute=0,second=0,microsecond=0)
date = date1.strftime("%Y%m%d")

if __name__  == "__main__":
    fileDir = rf"{Dir}\gfs{date}"
    if not os.path.exists(fileDir):
        os.mkdir(fileDir)
    if preci == 0.25:
        pre='0p25'
    elif preci == 0.5:
        pre='0p50'
    elif preci ==1:
        pre='1p00'
    gfs = f"https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_{pre}.pl?file=gfs.t00z.pgrb2.{pre}.f"
    for run in runs:
        url=f"{gfs}{run:03d}&all_lev=on&all_var=on&subregion=&leftlon={lon[0]}&rightlon={lon[1]}&toplat={lat[1]}&bottomlat={lat[0]}&dir=%2Fgfs.{date}%2F00%2Fatmos"
        file_dt = (date1+dt.timedelta(hours=int(run))).strftime("%Y%m%d_%H%M%S")
        print(url)
        filenm = wget.download(url, out=fileDir)
        print(filenm)