def rename_wrfout(path):
    import os
    FilePathDir=path
    for dirpath, dirnames, filenames in os.walk(FilePathDir):
        for file in filenames:
            if file.split('_')[0] == "wrfout":
                if file.split('.')[-1] != "nc":
                    filename, type = os.path.splitext(file)
                    os.rename(os.path.join(dirpath,file), os.path.join(dirpath, filename.replace('','') + '.nc') ) 

def wrf_plot(variable,datapath='',tm=[],filenm='',domain='d01',savepath='',suffix='',isave=False,ishow=True,var_lim='default',icontourf=True):
    # this function is used to plot variable field.
    # usage_example:
    # datapath=r"E:\nc_file\Data\INFA_Test\WRFRUN"
    # tm=[2021,7,25,12]
    # 'datapath' & 'tm' are here to determine the filename. You can also directely assign the filename by 'filenm'.
    # variable='wind10' or 'rain' or 'slp'
    # domain='d01' 'd02' etc.
    # savepath='wrf_result_analysis'
    # suffix='_infa_bogus'
    import matplotlib.pyplot as plt
    import numpy as np
    import cartopy.crs as ccrs   
    import cartopy.feature as cfeature
    from wrf import to_np,interplevel
    from matplotlib.pyplot import get_cmap
    from wrf import getvar, interpline, CoordPair,get_cartopy, cartopy_xlim, cartopy_ylim, latlon_coords,geo_bounds
    import netCDF4 as nc
    import matplotlib.ticker as ticker
    import matplotlib
    import os
    import datetime
    
    if filenm=='':
        wrffile=nc.Dataset(rf'{datapath}/wrfout_{domain}_{tm[0]:04d}-{tm[1]:02d}-{tm[2]:02d}_{tm[3]:02d}0000.nc')
    else:
        tmstr=filenm[-20:-3]
        date = datetime.datetime.strptime(tmstr, "%Y-%m-%d_%H%M%S")
        tm=date.timetuple()[:4]
        wrffile=nc.Dataset(filenm)
        
    U=getvar(wrffile,'U10')
    lats, lons = latlon_coords(U)
    cart_proj = get_cartopy(U)

    xlim=cartopy_xlim(U)
    ylim=cartopy_ylim(U)

    X, Y = np.meshgrid(to_np(lons), to_np(lats))
    if not savepath=='':
        savepath=savepath+'/'
    #################################
    if variable=='wind10':
        V=getvar(wrffile, 'V10')
        WIND=(U**2+V**2)**0.5
        
        fig1 = plt.figure(figsize=(8,6),dpi=200)
        ax1 = plt.axes(projection=cart_proj)
        ax1.add_feature(cfeature.LAND)
        ax1.coastlines('10m', linewidth=0.8)

        
        if not var_lim=='default':
            # norm = matplotlib.colors.Normalize(vmin=var_lim[0],vmax=var_lim[1])
            levels=list(np.linspace(var_lim[0],var_lim[1],num=10,dtype='int32',endpoint=False))
        else:
            levels=None
            
        if icontourf:
            plt.contourf(to_np(lons), to_np(lats), to_np(WIND), 10,
                         transform=ccrs.PlateCarree(),
                         cmap=drow_wind_heatmap,levels=levels,extend='both')#,norm=norm
            plt.colorbar(ax=ax1, shrink=.98)
        else:
            cs = plt.contour(to_np(lons), to_np(lats), to_np(WIND),10, colors="black",
                        transform=ccrs.PlateCarree(),levels=levels)
            ax1.clabel(cs,cs.levels, fontsize=8, colors='k')
        
        ax1.set_xlim(xlim)
        ax1.set_ylim(ylim)
        gl=ax1.gridlines(color="black", linestyle="dotted",draw_labels=True,dms=False,x_inline=False, y_inline=False,rotate_labels=False)
        gl.top_labels = False
        gl.right_lables = False
        plt.title(f"Wind Speed-10m (m/s)\ndatetime:{tm[0]:04d}-{tm[1]:02d}-{tm[2]:02d} {tm[3]:02d}:00:00")
        if isave:
            plt.savefig(f'{savepath}{tm[0]:04d}{tm[1]:02d}{tm[2]:02d}{tm[3]:02d}_Wind_Speed_10m{suffix}.png',dpi=200)
        if ishow:
            plt.show()
        else:
            plt.close()
        
        fig2 = plt.figure(figsize=(8,6),dpi=200)
        ax2 = plt.subplot(1,1,1,projection=cart_proj)
        ax2.add_feature(cfeature.LAND)
        ax2.coastlines('10m', linewidth=0.8)
        q = plt.quiver(to_np(lons), to_np(lats), to_np(U), to_np(V),transform=ccrs.PlateCarree())
        ax2.quiverkey(q, X=0.8, Y=0.9, U=20, label='20m/s ', labelpos='E',coordinates='figure')
        ax2.set_xlim(xlim)
        ax2.set_ylim(ylim)
        gl=ax2.gridlines(color="black", linestyle="dotted",draw_labels=True,dms=False,x_inline=False, y_inline=False,rotate_labels=False,alpha=0)
        gl.top_labels = False
        gl.right_lables = False
        plt.title(f"Wind Field\ndatetime:{tm[0]:04d}-{tm[1]:02d}-{tm[2]:02d} {tm[3]:02d}:00:00")
        if isave:
            plt.savefig(f'{savepath}{tm[0]:04d}{tm[1]:02d}{tm[2]:02d}{tm[3]:02d}_wind_filed_tc{suffix}.png',dpi=200)
        if ishow:
            plt.show()
        else:
            plt.close()
            
    elif variable=='slp':
        P=getvar(wrffile, 'slp')
        p=getvar(wrffile, 'pressure')
        z=getvar(wrffile, 'z')
        ht_500mb = interplevel(z, p, 500.)
        ht_500mb = ht_500mb/10.
        lon_500mb = ht_500mb.XLONG
        lat_500mb = ht_500mb.XLAT
        
        fig1 = plt.figure(figsize=(8,6),dpi=200)
        ax1 = plt.axes(projection=cart_proj)
        ax1.add_feature(cfeature.LAND)
        ax1.coastlines('10m', linewidth=0.8)
        if not var_lim=='default':
            levels=list(np.linspace(var_lim[0],var_lim[1],num=10,dtype='int32',endpoint=False))
        else:
            levels=None
        if icontourf:
            plt.contourf(to_np(lons), to_np(lats), to_np(P), 10,
                         transform=ccrs.PlateCarree(),
                         cmap=get_cmap("Spectral_r").reversed(),levels=levels,extend='both')#,norm=norm
            plt.colorbar(ax=ax1, shrink=.98)
        else:
            cs = plt.contour(to_np(lons), to_np(lats), to_np(P),10, colors="black",
                        transform=ccrs.PlateCarree(),levels=levels)
            ax1.clabel(cs,cs.levels, fontsize=8, colors='k')
        ax1.set_xlim(xlim)
        ax1.set_ylim(ylim)
        gl=ax1.gridlines(color="black", linestyle="dotted",draw_labels=True,dms=False,x_inline=False, y_inline=False,rotate_labels=False)
        gl.top_labels = False
        gl.right_lables = False
        plt.title(f"Sea Level Pressure (hPa)\ndatetime:{tm[0]:04d}-{tm[1]:02d}-{tm[2]:02d} {tm[3]:02d}:00:00")
        if isave:
            plt.savefig(f'{savepath}{tm[0]:04d}{tm[1]:02d}{tm[2]:02d}{tm[3]:02d}_Sea_Level_Pressure{suffix}.png',dpi=200)
        if ishow:
            plt.show()
        else:
            plt.close()
        
        fig2 = plt.figure(figsize=(8,6),dpi=200)
        ax2 = plt.axes(projection=cart_proj)
        ax2.add_feature(cfeature.LAND)
        ax2.coastlines('10m', linewidth=0.8)
        
        if not var_lim=='default':
            levels=list(np.linspace(var_lim[2],var_lim[3],num=10,dtype='int32',endpoint=False))
        else:
            levels=None
        if icontourf:
            plt.contourf(lon_500mb, lat_500mb, ht_500mb, 10,
                         transform=ccrs.PlateCarree(),
                         cmap=get_cmap("Spectral_r").reversed(),levels=levels,extend='both')
            plt.colorbar(ax=ax2, shrink=.98)
        else:
            cs = plt.contour(lon_500mb, lat_500mb, ht_500mb,10, colors="black",
                        transform=ccrs.PlateCarree(),levels=levels)
            ax2.clabel(cs,cs.levels, fontsize=8, colors='k')

        ax2.set_xlim(xlim)
        ax2.set_ylim(ylim)
        gl=ax2.gridlines(color="black", linestyle="dotted",draw_labels=True,dms=False,x_inline=False, y_inline=False,rotate_labels=False)
        gl.top_labels = False
        gl.right_lables = False
        plt.title(f"Geopotential height at 500hPa (dagpm)\ndatetime:{tm[0]:04d}-{tm[1]:02d}-{tm[2]:02d} {tm[3]:02d}:00:00")
        if isave:
            plt.savefig(f'{savepath}{tm[0]:04d}{tm[1]:02d}{tm[2]:02d}{tm[3]:02d}_Geopotential_height{suffix}.png',dpi=200)
        if ishow:
            plt.show()
        else:
            plt.close()
    elif variable=='rain':
        R1 = getvar(wrffile, 'RAINC')
        R2 = getvar(wrffile, 'RAINNC')
        R3 = getvar(wrffile, 'RAINSH')
        R=R1+R2+R3
        
        fig1 = plt.figure(figsize=(8,6),dpi=200)
        ax1 = plt.axes(projection=cart_proj)
        ax1.add_feature(cfeature.LAND)
        ax1.coastlines('10m', linewidth=0.8)
        if not var_lim=='default':
            levels=list(np.linspace(var_lim[0],var_lim[1],num=10,dtype='int32',endpoint=False))
        else:
            levels=None
        if icontourf:
            plt.contourf(to_np(lons), to_np(lats), to_np(R), 10,
                         transform=ccrs.PlateCarree(),
                         cmap=get_cmap("Spectral_r"),levels=levels,extend='both')#,norm=norm
            plt.colorbar(ax=ax1, shrink=.98)
        else:
            cs = plt.contour(to_np(lons), to_np(lats), to_np(R),10, colors="black",
                        transform=ccrs.PlateCarree(),levels=levels)
            ax1.clabel(cs,cs.levels, fontsize=8, colors='k')
        ax1.set_xlim(xlim)
        ax1.set_ylim(ylim)
        gl=ax1.gridlines(color="black", linestyle="dotted",draw_labels=True,dms=False,x_inline=False, y_inline=False,rotate_labels=False)
        gl.top_labels = False
        gl.right_lables = False
        plt.title(f"Rain (mm)\ndatetime:{tm[0]:04d}-{tm[1]:02d}-{tm[2]:02d} {tm[3]:02d}:00:00")
        if isave:
            plt.savefig(f'{savepath}{tm[0]:04d}{tm[1]:02d}{tm[2]:02d}{tm[3]:02d}_Rain{suffix}.png',dpi=200)
        if ishow:
            plt.show()
        else:
            plt.close()

def wrf_gif_plot(filepath,variable,domain='d01',start_index=0,end_index=0,time_interval_out=48,iremove=True,var_lim='default'):
    import glob
    import matplotlib.animation as anim
    import matplotlib.pyplot as plt
    import os
    import datetime
    import imageio

    fnms_all=glob.glob(filepath+'\wrfout_'+domain+'*.nc')
    times = ['' for _ in range(len(fnms_all))]
    i=0
    for fnm in fnms_all:
        times[i]=fnm[-20:-3]
        i+=1
    date0 = datetime.datetime.strptime(times[0], "%Y-%m-%d_%H%M%S")
    date1 = datetime.datetime.strptime(times[1], "%Y-%m-%d_%H%M%S")
    history_interval = (date1 - date0).total_seconds()/3600  #hours
    file_interval=round(time_interval_out/history_interval)
    if end_index==0:
        fnms=fnms_all[start_index::file_interval]
    else:
        fnms=fnms_all[start_index:end_index:file_interval]

    for fnm in fnms:
        wrf_plot(variable,filenm=fnm,ishow=False,isave=True,savepath='tmp',suffix='_tmp',var_lim=var_lim)
        print(fnm)

    def convert_img2gif(img_path,img_suffix,gif_name):
        files = glob.glob(f'{img_path}/*{img_suffix}')
        files.sort()
        frames = []
        for file in files:
            frames.append(imageio.imread(file))
        imageio.mimsave(gif_name, frames, 'GIF', duration = 500,loop=0)

    if variable=='wind10':
        convert_img2gif('tmp','Wind_Speed_10m_tmp.png','Wind_Speed_10m.gif')
        convert_img2gif('tmp','wind_filed_tc_tmp.png','wind_filed_tc.gif')
    elif variable=='slp':
        convert_img2gif('tmp','Sea_Level_Pressure_tmp.png','Sea_Level_Pressure.gif')
        convert_img2gif('tmp','Geopotential_height_tmp.png','Geopotential_height.gif')
    elif variable=='rain':
        convert_img2gif('tmp','Rain_tmp.png','Rain.gif')
    if iremove:
        tmpfiles=glob.glob('tmp/*tmp.png')
        for tmpfile in tmpfiles:
            os.remove(tmpfile)
    
def wrf_track(FilePathDir,lonTC0=-60.0,latTC0=-60.0,domain='d01',suffix='',ifbogus=False,track_file = 'data/cma_bst_infa.txt'):
    import numpy as np
    import netCDF4 as nc
    import datetime
    from wrf import getvar, ALL_TIMES, to_np
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    from cartopy.io.shapereader import Reader
    import cartopy.feature as cfeature
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
    import pandas as pd
    import xarray as xr
    import glob
    import os
    import openpyxl
    import pandas as pd
    
    def nearest_position( stn_lat, stn_lon, lat2d, lon2d):
        difflat = stn_lat - lat2d;
        difflon = stn_lon - lon2d;
        rad = np.multiply(difflat,difflat)+np.multiply(difflon , difflon)
        aa=np.where(rad==np.min(rad))
        ind=np.squeeze(np.array(aa))
        return tuple(ind)
    
    all_files=glob.glob(FilePathDir+f'\wrfout_{domain}*.nc')
    all_files.sort()
    
    if ifbogus:
        filenms=all_files[:]
    else:
        filenms=all_files[1:]
    times = ['' for _ in range(len(filenms))]
    i=0
    for filenm in filenms:
        times[i]=filenm[-20:-3]
        i+=1

    date0 = datetime.datetime.strptime(times[0], "%Y-%m-%d_%H%M%S")
    date1 = datetime.datetime.strptime(times[1], "%Y-%m-%d_%H%M%S")
    history_interval = (date1 - date0).total_seconds()/3600  #hours
    nt = len(times)
    for it in range(nt):
        time=times[it]
        wrfout      = nc.Dataset(filenms[it], mode="r")
        # initialize
        if it==0:
            slp    = to_np(getvar(wrfout ,"slp"  ,units="hPa" ))
            lat2D  = to_np(getvar(wrfout, "lat"  ))  # units: decimal degrees
            lon2D  = to_np(getvar(wrfout, "lon"  ))  # units: decimal degrees
            ny, nx = np.shape(lat2D)
            lonMax = np.max(lon2D)
            lonMin = np.min(lon2D)
            latMax = np.max(lat2D)
            latMin = np.min(lat2D)
            date = [] # 2020-07-20T00
            lons = [] # degree
            lats = [] # degree
            pmin = [] # hPa
            vmax = [] # m/s

            if latTC0 > -60.0 :
                latTC = latTC0
                lonTC = lonTC0
            else:
                slpMin = np.min(slp[:,:])
                indexMin = np.argwhere(slp[:,:] == slpMin)
                jTC = indexMin[0][0]
                iTC = indexMin[0][1]
                lonTC = lon2D[jTC,iTC]
                latTC = lat2D[jTC,iTC]

        slp         = to_np(getvar(wrfout ,"slp"         ,units="hPa"))
        wspd_wdir10 = to_np(getvar(wrfout ,"wspd_wdir10" ,units="m s-1"))
        wspd10      = wspd_wdir10[0]

        indexTC = nearest_position(latTC, lonTC, lat2D, lon2D)
        jTC = indexTC[0]
        iTC = indexTC[1]

        jTC = np.max((1,jTC))    # jTC [1,ny-2]
        jTC = np.min((jTC,ny-2))
        iTC = np.max((1,iTC))    # iTC [1,nx-2]
        iTC = np.min((iTC,nx-2))

        dLat = lat2D[jTC+1,iTC] - lat2D[jTC,iTC]
        dLon = lon2D[jTC,iTC+1] - lon2D[jTC,iTC]
        dAvg = (dLat + dLon)/2.0

        if abs(latTC) < 30.0 :
           radius = 0.5*history_interval  # 0.5 degree/hour
        else:
           radius = 1.0*history_interval  # 1.0 degree/hour
        if it==0 :
           radius = 0.5
        indexRadius = int(radius/dAvg) + 1

        iStart = iTC - indexRadius
        iEnd   = iTC + indexRadius
        jStart = jTC - indexRadius
        jEnd   = jTC + indexRadius
        jStart = np.max((1,jStart))
        jEnd   = np.min((jEnd,ny-2))
        iStart = np.max((1,iStart))
        iEnd   = np.min((iEnd,nx-2))

        slpMin = np.min(slp[jStart:jEnd,iStart:iEnd])
        w10Max = np.max(wspd10[jStart:jEnd,iStart:iEnd])
        indexMin = np.argwhere(slp[jStart:jEnd,iStart:iEnd] == slpMin)
        jTC = indexMin[0][0] + jStart
        iTC = indexMin[0][1] + iStart
        lonTC = lon2D[jTC,iTC]
        latTC = lat2D[jTC,iTC]
        print("date:", str(times[it])[0:19],"TC center:",round(lonTC,2), round(latTC,2)," p min:",round(slpMin,2), " vmax:",round(w10Max,2))
        date.append(str(times[it])[0:19])
        lons.append(round(lonTC,2))
        lats.append(round(latTC,2))
        pmin.append(round(slpMin,2))
        vmax.append(round(w10Max,2))
    
    # save tropycal center coordinates, pmin and vmax
    tm_dt=[datetime.datetime.strptime(date[i],"%Y-%m-%d_%H%M%S") for i in range(len(date))]
    tm_str2=[tm_dt[i].strftime('%Y%m%d%H') for i in range(len(tm_dt))]
    mat = {'time':tm_str2,'lon':lons,'lat':lats,'mslp':pmin,'vmax':vmax}
    mat_df=pd.DataFrame(mat)
    mat_df.to_csv('wrfout_track'+suffix+'.txt',index=False)
    
    #read CMA-STI best track
    f_Con = track_file
    col_names =['date','grade','lat', 'lon', 'pres', 'vmax']
    widths = [10,2,4,5,4,3]
    df = pd.read_fwf(f_Con,usecols=[0,1,2,3,4,5],widths=widths,names=col_names)
    latObs  = df['lat'].values[:]
    lonObs  = df['lon'].values[:]
    latObs  = np.array(latObs)/10
    lonObs  = np.array(lonObs)/10
    
    ### plot track
    fig = plt.figure()
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([lonMin, lonMax, latMin, latMax],crs=ccrs.PlateCarree())

    ax.add_feature(cfeature.LAND)
    gl = ax.gridlines(color="black", linestyle="dotted",draw_labels=True,dms=False,x_inline=False, y_inline=False,rotate_labels=False)
    gl.top_labels = False
    gl.right_lables = False
    ax.coastlines()
    for it in range(1,len(lons)):
        if it == 1:
           plt.plot((lons[it-1],lons[it]), (lats[it-1],lats[it]),color='red',linewidth=1.5,transform=ccrs.PlateCarree(), label="wrf")
        else:
           plt.plot((lons[it-1],lons[it]), (lats[it-1],lats[it]),color='red',linewidth=1.5,transform=ccrs.PlateCarree())
    for it in range(1,len(lonObs)):
        if it == 1:
            plt.plot((lonObs[it-1],lonObs[it]), (latObs[it-1],latObs[it]),color='black',linewidth=2,transform=ccrs.PlateCarree(), label="obs")
        else:
            plt.plot((lonObs[it-1],lonObs[it]), (latObs[it-1],latObs[it]),color='black',linewidth=2,transform=ccrs.PlateCarree())
    plt.legend()
    plt.savefig('typhoon tracks'+suffix+'.png',dpi=200)
    plt.show()
    
    
def fnl_plot(tm,wrfpath,fnlpath,domain='d01'):
    # you need to give a wrfpath (with wrfout files renamed by function rename_wrfout) to determine the extent of target fnl plot
    import matplotlib.pyplot as plt
    import numpy as np
    import cartopy.crs as ccrs   
    import cartopy.feature as cfeature
    from wrf import to_np,interplevel
    from matplotlib.pyplot import get_cmap
    from wrf import getvar, interpline, CoordPair,get_cartopy, cartopy_xlim, cartopy_ylim, latlon_coords,geo_bounds
    import netCDF4 as nc
    import matplotlib.ticker as ticker
    import os
    import datetime
    import cfgrib

    wrffile=nc.Dataset(rf'{wrfpath}/wrfout_{domain}_{tm[0]:04d}-{tm[1]:02d}-{tm[2]:02d}_{tm[3]:02d}0000.nc')
    U=getvar(wrffile,'U10')
    lats, lons = latlon_coords(U)
    cart_proj = get_cartopy(U)
    xlim=cartopy_xlim(U)
    ylim=cartopy_ylim(U)
    
    ds_fnl=cfgrib.open_datasets(rf"{fnlpath}/fnl_{tm[0]:04d}{tm[1]:02d}{tm[2]:02d}_{tm[3]:02d}_00.grib2")
    
    U_fnl=ds_fnl[7]['u10'].values
    V_fnl=ds_fnl[7]['v10'].values
    WIND_fnl=(U_fnl**2+V_fnl**2)**0.5
    lons_fnl=ds_fnl[7].longitude.values
    lats_fnl=ds_fnl[7].latitude.values
    X_fnl, Y_fnl = np.meshgrid(to_np(lons_fnl), to_np(lats_fnl))
    fig = plt.figure(figsize=(8,6),dpi=200)
    ax = plt.subplot(1,1,1,projection=cart_proj)
    ax.add_feature(cfeature.LAND)
    ax.coastlines('10m', linewidth=0.8)
    q = plt.quiver(to_np(lons_fnl), to_np(lats_fnl), to_np(U_fnl), to_np(V_fnl),transform=ccrs.PlateCarree())
    ax.quiverkey(q, X=0.7, Y=0.9, U=20, label='20m/s ', labelpos='E',coordinates='figure')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    gl=ax.gridlines(color="black", linestyle="dotted",draw_labels=True,dms=False,x_inline=False, y_inline=False,rotate_labels=False,alpha=0)
    gl.top_labels = False
    gl.right_lables = False
    plt.title(f"Wind Field  datetime:{tm[0]:04d}-{tm[1]:02d}-{tm[2]:02d} {tm[3]:02d}:00:00")
    plt.savefig(f'{tm[0]:04d}{tm[1]:02d}{tm[2]:02d}{tm[3]:02d}_wind_filed_original_fnl.png',dpi=200)
    plt.show()
    
def difference_between_wrfout(variable,paths,tm,plot_title='',domain='d01',isave=True,ishow=True):
    '''
        usage example:
        import sys
        sys.path.append(r"D:\Anaconda\myPythonProj")
        import WRFpost

        path1=r"E:\nc_file\Data\INFA_Test\WRFRUN"
        path2=r"E:\nc_file\Data\INFA_Test\WRFRUN_bogus"
        tm=[2021,7,25,12]
        WRFpost.difference_between_wrfout('wind10',[path1,path2],tm)
    '''
    import matplotlib.pyplot as plt
    import numpy as np
    import cartopy.crs as ccrs   
    import cartopy.feature as cfeature
    from matplotlib.pyplot import get_cmap
    from wrf import to_np,getvar,get_cartopy, cartopy_xlim, cartopy_ylim, latlon_coords
    import netCDF4 as nc
    
    # 2 wrfout files must have same grids/meshes
    fnm1=rf'{paths[0]}/wrfout_{domain}_{tm[0]:04d}-{tm[1]:02d}-{tm[2]:02d}_{tm[3]:02d}0000.nc'
    fnm2=rf'{paths[1]}/wrfout_{domain}_{tm[0]:04d}-{tm[1]:02d}-{tm[2]:02d}_{tm[3]:02d}0000.nc'
    wrffile1=nc.Dataset(fnm1)
    wrffile2=nc.Dataset(fnm2)

    U=getvar(wrffile1,'U10')
    lats, lons = latlon_coords(U)
    cart_proj = get_cartopy(U)
    xlim=cartopy_xlim(U)
    ylim=cartopy_ylim(U)
    X, Y = np.meshgrid(to_np(lons), to_np(lats))
    if variable=='wind10':
        U1 = getvar(wrffile1, 'U10')
        U2 = getvar(wrffile2, 'U10')
        V1 = getvar(wrffile1, 'V10')
        V2 = getvar(wrffile2, 'V10')
        VAR1=(U1**2+V1**2)**0.5
        VAR2=(U2**2+V2**2)**0.5
    elif variable=='slp':
        VAR1 = getvar(wrffile1, 'slp')
        VAR2 = getvar(wrffile2, 'slp')
    elif variable=='rain':
        R1 = getvar(wrffile1, 'RAINC')
        R2 = getvar(wrffile1, 'RAINNC')
        R3 = getvar(wrffile1, 'RAINSH')
        VAR1=R1+R2+R3
        R1 = getvar(wrffile2, 'RAINC')
        R2 = getvar(wrffile2, 'RAINNC')
        R3 = getvar(wrffile2, 'RAINSH')
        VAR2=R1+R2+R3

    fig1 = plt.figure(figsize=(8,6),dpi=200)
    ax1 = plt.axes(projection=cart_proj)
    ax1.add_feature(cfeature.LAND)
    ax1.coastlines('10m', linewidth=0.8)


    plt.contourf(to_np(lons), to_np(lats), to_np(VAR2-VAR1), 10,
                 transform=ccrs.PlateCarree(),
                 cmap=get_cmap("Spectral_r"),extend='both')#,norm=norm
    plt.colorbar(ax=ax1, shrink=.98)

    ax1.set_xlim(xlim)
    ax1.set_ylim(ylim)
    gl=ax1.gridlines(color="black", linestyle="dotted",draw_labels=True,dms=False,x_inline=False, y_inline=False,rotate_labels=False)
    gl.top_labels = False
    gl.right_lables = False
    if plot_title=='':
        plot_title=f"{variable} difference (m/s)\ndatetime:{tm[0]:04d}-{tm[1]:02d}-{tm[2]:02d} {tm[3]:02d}:00:00"
    plt.title(plot_title)
    if isave:
        plt.savefig(f'{variable} difference.png',dpi=200)
    if ishow:
        plt.show()
    else:
        plt.close()
        
def difference_track(wrf_track_fnm,wrfoutnm,storm_name,storm_year,ibtracs_fnm=r'E:\OneDrive\OneDrive - stu.ouc.edu.cn\document\12_软件学习\WRF\ibtracs.ALL.list.v04r00.csv'):
    '''
        usage example:
        wrf_track_fnm='wrfout_track.txt'
        wrfoutnm=r"E:\nc_file\Data\INFA_Test\WRFRUN\wrfout_d01_2021-07-22_120000.nc" #给出一个wrfout文件用于确定坐标系以及绘制地图的范围
        WRFpost.difference_track(wrf_track_fnm,wrfoutnm,'IN-FA',2021)
    '''
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    from tropycal import tracks, utils
    import numpy as np
    import datetime as dt
    import sys
    sys.path.append(r"D:\Anaconda\myPythonProj")
    import my_utils
    import pandas as pd

    ibtracs = tracks.TrackDataset(ibtracs_url=ibtracs_fnm,basin='all',source='ibtracs',ibtracs_mode='jtwc_neumann',catarina=True)
    storm_ibt = ibtracs.get_storm((storm_name,storm_year)).interp()
    tm_ibt=storm_ibt['time']
    lat_ibt=storm_ibt['lat']
    lon_ibt=storm_ibt['lon']
    vmax_ibt=storm_ibt['vmax']*0.5144
    slp_ibt=storm_ibt['mslp']

    storm_dict=utils.create_storm_dict(wrf_track_fnm,storm_name+'_simu','001')
    storm_sim=tracks.Storm(storm_dict).interp()
    tm_sim=storm_sim['time']
    lat_sim=storm_sim['lat']
    lon_sim=storm_sim['lon']
    vmax_sim=storm_sim['vmax']
    slp_sim=storm_sim['mslp']

    my_colors=my_utils.my_colors('contrast',n_colors=10)
    fig=plt.subplots(1,3,figsize=(18,6),dpi=600,tight_layout=True)#,subplot_kw={'projection':ccrs.PlateCarree()}

    # ax1: track plot
    from wrf import getvar, ALL_TIMES, to_np
    import netCDF4 as nc
    wrfout = nc.Dataset(wrfoutnm, mode="r")
    lat2D  = to_np(getvar(wrfout, "lat"  ))  # units: decimal degrees
    lon2D  = to_np(getvar(wrfout, "lon"  ))  # units: decimal degrees
    lonMax = np.max(lon2D)
    lonMin = np.min(lon2D)
    latMax = np.max(lat2D)
    latMin = np.min(lat2D)
    ax1=plt.subplot(131,projection=ccrs.PlateCarree())
    ax1.set_extent([lonMin, lonMax, latMin, latMax],crs=ccrs.PlateCarree())
    ax1.add_feature(cfeature.LAND)
    ax1.plot(lon_ibt,lat_ibt,color=my_colors[0],linewidth=1.5,linestyle='--',transform=ccrs.PlateCarree(), label="ibtracs")
    ax1.plot(lon_sim,lat_sim,color=my_colors[1],linewidth=1.5,transform=ccrs.PlateCarree(), label="wrf")
    gl = ax1.gridlines(color="black", linestyle="dotted",draw_labels=True,dms=False,x_inline=False, y_inline=False,rotate_labels=False)
    gl.top_labels = False
    gl.right_lables = False
    ax1.coastlines()
    ax1.legend()
    # ax2: vmax time series plot
    ax2=plt.subplot(132)
    plt.xticks(rotation=45)
    ax2.plot(tm_ibt,vmax_ibt,color=my_colors[0],linewidth=1.5,linestyle='--', label="ibtracs")
    ax2.plot(tm_sim,vmax_sim,color=my_colors[1],linewidth=1.5, label="wrf")
    ax2.set_ylabel('vmax m/s')
    ax2.legend()
    # ax3: pmin time series plot
    ax3=plt.subplot(133)
    plt.xticks(rotation=45)
    ax3.plot(tm_ibt,slp_ibt,color=my_colors[0],linewidth=1.5,linestyle='--', label="ibtracs")
    ax3.plot(tm_sim,slp_sim,color=my_colors[1],linewidth=1.5, label="wrf")
    ax3.set_ylabel('pmin hPa')
    ax3.legend()
    plt.savefig(f'typhoon_tracks_{storm_name}.png',dpi=600)
    dic={"time":tm_ibt,"lon":lon_ibt,"lat":lat_ibt,"vmax":vmax_ibt,"pmin":slp_ibt}
    ibt_tropical=pd.DataFrame(dic)
    return ibtracs,ibt_tropical

def ts_plot(datapath,location=[],measure_fnm='',fnlpath='',domain='d01',suffix='',isave=True):
    '''
    this function is used to plot time series.
    usage_example:
        import datetime
        start=datetime.datetime.now()
        datapath=r"E:\nc_file\Data\INFA_Test\WRFRUN"
        location=[121.96,31.51]#崇明东滩湿地气象站？default: location=[]
        measure_fnm=r'E:\OneDrive\OneDrive - stu.ouc.edu.cn\document\13_我的野外数据\崇明气象数据\气象局20210501-20210823.xls'
        domain='d02'
        fnlpath=r'E:\nc_file\Data\INFA_Test'
        
        %load_ext autoreload
        %autoreload 2
        import sys
        sys.path.append(r"D:\Anaconda\myPythonProj")
        import WRFpost
        WRFpost.ts_plot(datapath,location,measure_fnm=measure_fnm,fnlpath=r'E:\nc_file\Data\INFA_Test',domain=domain)
        end=datetime.datetime.now()
        interval_time=(end-start).total_seconds()
        print(f'running time: {interval_time:6f}')
    '''
    
    import sys
    sys.path.append(r"D:\Anaconda\myPythonProj")
    import my_utils
    from scipy.interpolate import griddata
    import matplotlib.pyplot as plt
    import numpy as np
    import matplotlib
    import os
    import datetime
    import glob
    from salem.utils import get_demo_file
    import salem
    from matplotlib.lines import Line2D
    from wrf import to_np,getvar

    all_files=glob.glob(datapath+f'\wrfout_{domain}*.nc')
    all_files.sort()
    filenms=all_files[:]
    times = ['' for _ in range(len(filenms))]
    i=0
    for filenm in filenms:
        times[i]=filenm[-20:-3]
        i+=1
    datetimes = [datetime.datetime.now for _ in range(len(filenms))]
    for i,time in enumerate(times):
        datetimes[i]=datetime.datetime.strptime(time, "%Y-%m-%d_%H%M%S")

    date0 = datetimes[0]
    date1 = datetimes[1]
    history_interval = (date1 - date0).total_seconds()/3600  #hours
    nt = len(datetimes)

    for it in range(nt):
        time=datetimes[it]
        wrfout      = salem.open_wrf_dataset(filenms[it])
        lat2D  = wrfout.lat.values  # units: decimal degrees
        lon2D  = wrfout.lon.values  # units: decimal degrees
        ny, nx = np.shape(lat2D)
        u         = np.squeeze(wrfout.U10.values)
        v         = np.squeeze(wrfout.V10.values)
        t2        = np.squeeze(wrfout.T2.values)
        wind_speed = (u ** 2 + v ** 2) ** 0.5
        wind_direction = np.mod(180 + (180 / np.pi) * np.arctan2(v, u),360) # units: degrees
        t2=t2-273.15 #units: degree celsius
        if it==0:
            # initialize
            ws = [] # m/s
            wd = [] # degree
            temp = [] # degree Celsius
            ws3D=np.zeros((1,ny,nx))
            wd3D=np.zeros((1,ny,nx))
            t23D=np.zeros((1,ny,nx))
            
        if not location==[]:
            lon1D=np.squeeze(lon2D.reshape(1,-1))
            lat1D=np.squeeze(lat2D.reshape(1,-1))
            wind_speed1D=np.squeeze(wind_speed.reshape(1,-1))
            wind_direction1D=np.squeeze(wind_direction.reshape(1,-1))
            t21D=np.squeeze(t2.reshape(1,-1))
            # interpolate
            wsi = griddata((lon1D,lat1D), wind_speed1D, (location[0], location[1]), method='linear').tolist()
            wdi = griddata((lon1D,lat1D), wind_direction1D, (location[0], location[1]), method='linear').tolist()
            t2i = griddata((lon1D,lat1D), t21D, (location[0], location[1]), method='linear').tolist()
            # append
            ws.append(wsi)
            wd.append(wdi)
            temp.append(t2i)
        else:
            ws3D=np.concatenate((ws3D,wind_speed[np.newaxis, :,:]), axis=0)
            wd3D=np.concatenate((wd3D,wind_direction[np.newaxis, :,:]), axis=0)
            t23D=np.concatenate((t23D,t2[np.newaxis, :,:]), axis=0)
            
    if not fnlpath=='':
        import cfgrib
        all_fnlnm=glob.glob(fnlpath+rf'\fnl_*.grib2')
        all_fnlnm.sort()
        fnlnms=all_fnlnm[:]
        fnltimes = ['' for _ in range(len(fnlnms))]
        for i,fnlnm in enumerate(fnlnms):
            fnltimes[i]=fnlnm[-20:-6]
        fnldts = [datetime.datetime.now for _ in range(len(fnlnms))]
        for i,time in enumerate(fnltimes):
            fnldts[i]=datetime.datetime.strptime(time, "%Y%m%d_%H_%M")

        for it,fnlnm in enumerate(fnlnms):
            ds_fnl=cfgrib.open_datasets(fnlnm)
            u=ds_fnl[7]['u10'].values
            v=ds_fnl[7]['v10'].values
            t2=ds_fnl[8]['t2m'].values
            wind_speed=(u**2+v**2)**0.5
            wind_direction = np.mod(180 + (180 / np.pi) * np.arctan2(v, u),360) # units: degrees
            t2=t2-273.15 #units: degree celsius

            lons_fnl=ds_fnl[7].longitude.values
            lats_fnl=ds_fnl[7].latitude.values
            lon2D,lat2D = np.meshgrid(to_np(lons_fnl), to_np(lats_fnl))
            ny,nx=np.shape(lat2D)
            if it==0:
                # initialize
                fnl_ws = [] # m/s
                fnl_wd = [] # degree
                fnl_temp = [] # degree Celsius

            if not location==[]:
                lon1D=np.squeeze(lon2D.reshape(1,-1))
                lat1D=np.squeeze(lat2D.reshape(1,-1))
                wind_speed1D=np.squeeze(wind_speed.reshape(1,-1))
                wind_direction1D=np.squeeze(wind_direction.reshape(1,-1))
                t21D=np.squeeze(t2.reshape(1,-1))
                # interpolate
                wsi = griddata((lon1D,lat1D), wind_speed1D, (location[0], location[1]), method='linear').tolist()
                wdi = griddata((lon1D,lat1D), wind_direction1D, (location[0], location[1]), method='linear').tolist()
                t2i = griddata((lon1D,lat1D), t21D, (location[0], location[1]), method='linear').tolist()
                # append
                fnl_ws.append(wsi)
                fnl_wd.append(wdi)
                fnl_temp.append(t2i)

    my_colors=my_utils.my_colors('contrast',n_colors=10)
    if not location==[]:
        ws=np.array(ws)#*0.82-0.07
        wd=np.array(wd)
        uu=-ws*np.sin(wd*np.pi/180)
        vv=-ws*np.cos(wd*np.pi/180)
        
        import pandas as pd
        measurement=pd.read_excel(measure_fnm)
        ob_tm_str=measurement['观测时间'].values
        ob_tm_all = [datetime.datetime.now for _ in range(len(ob_tm_str))]
        for i,time in enumerate(ob_tm_str):
            ob_tm_all[i]=datetime.datetime.strptime(time, "%Y-%m-%d %H")-datetime.timedelta(hours=8)# convert Beijing Time
        # ob_tm_all=measurement['Time'].values.astype('M8[ms]').astype('O')
        ob_ind=my_utils.in_range(ob_tm_all,min(datetimes),max(datetimes))
        ob_tm=np.array(ob_tm_all)[ob_ind]
        ob_ws=np.array(measurement['2分钟平均风速'].values[ob_ind],dtype='float')
        ob_wd=np.array(measurement['2分钟平均风向'].values[ob_ind],dtype='float')
        ob_t2=np.array(measurement['气温'].values[ob_ind],dtype='float')
        # ob_ws=np.array(measurement['Wind Speed'].values[ob_ind],dtype='float')
        # ob_wd=np.array(measurement['Wind Direction'].values[ob_ind],dtype='float')
        # ob_t2=np.array(measurement['Temperature'].values[ob_ind],dtype='float')
        ob_uu=-ob_ws*np.sin(ob_wd*np.pi/180)
        ob_vv=-ob_ws*np.cos(ob_wd*np.pi/180)   
        
        fig = plt.figure(dpi=200, figsize=(10,6))
        fig.autofmt_xdate()
        plt.style.use('seaborn-v0_8-ticks')
        q1=plt.quiver(datetimes,np.zeros_like(ws),uu,vv,color= my_colors[6],width=0.0025,scale=120)#,scale=1,scale_units='height',angles='xy',width=0.0025
        q1.set_linestyle(':')
        q2=plt.quiver(ob_tm,np.zeros_like(ob_ws),ob_uu,ob_vv,color= my_colors[0],width=0.0025,scale=120)#
        plt.quiverkey(q1,0.15, 0.75, 20, '20m/s\n\nwrf', labelpos='N')
        plt.quiverkey(q2,0.15, 0.65, 20, 'obs', labelpos='N')
        if not fnlpath==[]:
            fnl_ws=np.array(fnl_ws)#*0.82-0.07
            fnl_wd=np.array(fnl_wd)
            fnl_temp=np.array(fnl_temp)
            fnl_uu=-fnl_ws*np.sin(fnl_wd*np.pi/180)
            fnl_vv=-fnl_ws*np.cos(fnl_wd*np.pi/180)
            q3=plt.quiver(fnldts,np.zeros_like(fnl_ws),fnl_uu,fnl_vv,color=my_colors[3],width=0.0025,scale=120)
            plt.quiverkey(q3,0.15, 0.55, 20, 'fnl', labelpos='N')
        if isave:
            plt.savefig(f'ts_wind{suffix}.png',dpi=600)
            import pandas as pd
            dic1={"wrf_tm":datetimes,"wrf_ws":ws,"wrf_wd":wd,"wrf_t2":temp}
            dic2={"ob_tm":ob_tm,"obs_ws":ob_ws,"obs_wd":ob_wd,"obs_t2":ob_t2}
            df1=pd.DataFrame(dic1)
            df2=pd.DataFrame(dic2)
            if not fnlpath==[]:
                dic3={"fnl_tm":fnldts,"fnl_ws":fnl_ws,"fnl_wd":fnl_wd,"fnl_t2":fnl_temp}
                df3=pd.DataFrame(dic3)
            with pd.ExcelWriter(f'time_series{suffix}.xlsx') as writer:
                df1.to_excel(writer,sheet_name='Sheet1',index=False)
                df2.to_excel(writer,sheet_name='Sheet2',index=False)
                if not fnlpath==[]:
                    df3.to_excel(writer,sheet_name='Sheet3',index=False)
    else:
        # Plot wind speed for each grid point as a line
        fig, (ax1, ax2, ax3) = plt.subplots(1,3,figsize=(18,6),dpi=200)
        fig.autofmt_xdate()
        for i in range(ws3D.shape[1]):
            for j in range(ws3D.shape[2]):
                ax1.plot(datetimes, ws3D[1:, i, j], color='Grey', linewidth=0.1)
                ax2.plot(datetimes, wd3D[1:, i, j], color='Grey', linewidth=0.1)
                ax3.plot(datetimes, t23D[1:, i, j], color='Grey', linewidth=0.1)

        # Calculate and plot the mean wind speed across all grid points for each time step
        ws_mean = ws3D.mean(axis=(1, 2))
        wd_mean = wd3D.mean(axis=(1, 2))
        t2_mean = t23D.mean(axis=(1, 2))
        
        lines1 = ax1.plot(datetimes, ws_mean[1:], linestyle='--', color=my_colors[0], linewidth=2)
        lines2 = ax2.plot(datetimes, wd_mean[1:], linestyle='--', color=my_colors[0], linewidth=2)
        lines3 = ax3.plot(datetimes, t2_mean[1:], linestyle='--', color=my_colors[0], linewidth=2)

        # Set x-axis label
        ax1.set_xlabel('Time')
        ax2.set_xlabel('Time')
        ax3.set_xlabel('Time')

        # Set y-axis label
        ax1.set_ylabel('Wind Speed (m/s)')
        ax2.set_ylabel('Wind Direction (degree)')
        ax3.set_ylabel('Temperature (degree calsius)')

        # Add legend
        legend_elements1 = [Line2D([0], [0], color='gray', lw=0.5, label='Wind speed for each grid point'),
                           Line2D([0], [0], linestyle='--', color=my_colors[0], lw=2, label='Mean wind speed across all grid points')]
        ax1.legend(handles=legend_elements1, loc='upper left')
        legend_elements2 = [Line2D([0], [0], color='gray', lw=0.5, label='Wind direction for each grid point'),
                   Line2D([0], [0], linestyle='--', color=my_colors[0], lw=2, label='Mean wind direction across all grid points')]
        ax2.legend(handles=legend_elements2, loc='upper left')
        legend_elements3 = [Line2D([0], [0], color='gray', lw=0.5, label='Temperature for each grid point'),
                   Line2D([0], [0], linestyle='--', color=my_colors[0], lw=2, label='Mean temperature across all grid points')]
        ax3.legend(handles=legend_elements3, loc='upper left')
        if isave:
            plt.savefig(f'ts_all_grid{suffix}.png',dpi=600)
        plt.show()
