def adv_post1(adv_array=None,dat_file=r"E:\OBchongming202306-09\ADV\N5783\CM062010.dat",
              start_time=None,ene_time=None,
              SNR_thr=5, cor_thr=70,save_fname=None,ireturn=False):
    '''
        这个函数用于将adv数据（3个方向的流速、压力）处理成burst平均，数据质量根据信噪比和相关系数进行筛选。如果一个burst中没有可用的数据则输出nan值。
        必要参数：dat_file
        可选参数：adv_array，如已通过np.loadtxt读取矩阵，则可通过这个参数避免重复读取
                  start_time，长度为6的tuple，如（2016,4,12,0,0,0）
                  end_time同上
                  SNR_thr和cor_thr表示区分数据质量的阈值，越高数据质量越好，但能筛选出的数据也就越少。
                  save_fname可以指定处理后的数据存储的文件名，ascii格式，默认以后缀_post.dat存储在原adv文件夹中。
                  ireturn是指是否返回结果矩阵。
    '''
    import numpy as np
    import datetime

    # initialize
    if adv_array==None:
        adv_array=np.loadtxt(dat_file)
    if save_fname==None:
        save_fname=dat_file[:-4]+"_post1.dat"
    hdr_file=dat_file[:-3]+"hdr"
    with open(hdr_file,'r') as f:
        lines=f.readlines()
        ststr=lines[6].split()
        etstr=lines[7].split()
        intstr=lines[13].split()
        spbstr=lines[14].split()
    burst_interval=float(intstr[2])
    samples_per_burst=int(spbstr[-1])
    if start_time==None:
        st=datetime.datetime.strptime(ststr[4]+" "+ststr[5],'%Y/%m/%d %H:%M:%S')
        et=datetime.datetime.strptime(etstr[4]+" "+etstr[5],'%Y/%m/%d %H:%M:%S')
    else:
        st=datetime.datetime(*start_time)
        et=datetime.datetime(*end_time)

    dt=datetime.timedelta(seconds=burst_interval)
    ts=np.arange(st,et,dt)
    if len(ts)==adv_array.shape[0]//samples_per_burst:
        len_ts=len(ts)
    elif len(ts)>adv_array.shape[0]//samples_per_burst:
        len_ts=adv_array.shape[0]//samples_per_burst
        ts=ts[:len_ts]
    else:
        print("请检查时间和采样频率参数")
        return
        
    # processing
    nan_index=[[],[],[]]
    for i in range(len_ts):
        sign=np.zeros(3)
        for j in range(samples_per_burst):
            index=i*samples_per_burst+j    
            for k in [8,9,10]:
                if adv_array[index,k]<SNR_thr or adv_array[index,k+3]<cor_thr:
                    adv_array[index,k-6]=np.nan # E-vel N-vel Up-vel
                    adv_array[index,14]=np.nan # pressure(dbar)
                    if adv_array[index,k]<SNR_thr:
                        sign[k-8]+=1
        for m in range(3):
            if (sign[m]/samples_per_burst)>0.2:
                nan_index[m].append(i)
                    
    new_mat=np.zeros((len_ts,4)) # E-vel N-vel Up-vel pressure(dbar)

    for i in range(len_ts):
        index=np.arange(i*samples_per_burst,(i+1)*samples_per_burst)
        for j,k in enumerate([2,3,4,14]):
            new_mat[i,j]=np.nanmean(adv_array[index,k])
    for m in range(3):
        new_mat[nan_index[m],m]=np.nan
    np.savetxt(save_fname,new_mat,fmt='% 1.7f',delimiter='    ')
    if ireturn:
        return new_mat,ts
        
        
def wave_post1(dataframe=None,wave_fnm=r"E:\OBchongming202306-09\浪潮仪wave\207175_20230831_0016.xlsx",
               freq=16,save_fname1=None,save_fname=None,iplot=False,ireturn=False):
    '''
        这个函数用于将水深数据成每秒平均
        必要参数：wave_fnm
        可选参数：dataframe，如已通过pd.DataFrame(columns=['Time', 'Depth'])读取矩阵，则可通过这个参数避免重复读取
                  save_fname和save_fname1可以指定处理后的数据存储的文件名，ascii格式，默认以后缀_depth_all.dat和_depth_per_sec.dat存储在原文件夹中。
                  iplot表示是否绘制水深时间序列图
                  ireturn是指是否返回结果矩阵。
    '''
    import numpy as np
    import datetime
    import pandas as pd
    
    # initialize
    if dataframe==None:
        xls = pd.ExcelFile(wave_fnm)
        dataframe= pd.DataFrame(columns=['Time', 'Depth'])

        # 遍历每个sheet并进行处理
        for sheet_name in xls.sheet_names[3:]:
            tem=pd.read_excel(wave_fnm,sheet_name=sheet_name,header=1).drop(["Pressure","Sea pressure"],axis=1)
            dataframe=dataframe.append(tem, ignore_index=True)
            print(f"Data from sheet '{sheet_name}' has been appended.")
    if save_fname==None:
        save_fname=wave_fnm[:-5]+'_depth_per_sec.dat'
        save_fname1=wave_fnm[:-5]+'_depth_all.dat'

    dep1=dataframe["Depth"].values[0:dataframe.shape[0]//freq*freq]
    tm=dataframe["Time"][::freq][:-1]
    
    np.savetxt(save_fname1,dep1,fmt='% 1.7f',delimiter='    ')
    print(f"{save_fname1}\n起始时间:{tm.values[0]}，间隔为{1/freq}秒")

    mean_dep=np.mean(dep1.reshape(-1,freq),axis=1)
    if iplot:
        from matplotlib import pyplot as plt
        fig,ax=plt.subplots(figsize=(18,3),tight_layout=True)
        ax.plot(tm,mean_dep)

    np.savetxt(save_fname,mean_dep,fmt='% 1.7f',delimiter='    ')
    print(f"{save_fname}\n起始时间:{tm.values[0]}，结束时间：{tm.values[-1]}，间隔为1秒")
    if ireturn:
        return dep1,mean_dep,tm
        
# 截取adv或wave数据的一部分：
def selecte_data(times_of_data,data,start_time=(2023,8,10,0,0,0),end_time=(2023,10,13,0,0,0)):
    # times_of_data格式为np.array(dtype=datetime.datetime)
    st=datetime.datetime(*start_time)
    et=datetime.datetime(*end_time)
    result=data[np.logical_and(times_of_data>=st, times_of_data<=et)]
    return result
    
def wave_post2(dat_fname,Fs=16,data_array=None,wave_f_max=8,wave_f_min=0.001,filter_n=12,time_group=600,save_fname1=None,save_fname2=None,ireturn=False): 
    import matplotlib.pyplot as plt
    import numpy as np
    import numpy.fft as fft
    import os

    if data_array is None:
        data=np.loadtxt(dat_fname)
    else:
        data=data_array.copy()

    T=1/Fs #采样周期
    L=len(data) #信号长度
    t=np.array([i*T for i in range(L)])
    data1=data.copy()
    for i in range(filter_n):
        complex_array=fft.fft(data1)
        freqs = fft.fftfreq(t.size, t[1] -t[0])#复数的模为信号的振幅（能量大小）

        filter_complex_wave=complex_array.copy()
        filter_complex_wave[np.logical_or(freqs>wave_f_max,freqs<wave_f_min)]=0
        data1=fft.ifft(filter_complex_wave).real#S_new是ifft变换后的序列
    
    S_ifft_wave=data1
    dep=data-S_ifft_wave
    S_ifft_wave[dep<0]=np.nan
    dep[dep<0]=np.nan

    n=S_ifft_wave.shape[0]%(Fs*time_group) #规定每10min（600s）为一组，去掉最后面多余的数据

    if n==0:
        n=None
    S_ifft_wave=S_ifft_wave[:-n]
    
    if save_fname1 is None:
        fdir,_=os.path.split(dat_fname)
        save_fname1=fdir+'/'+'wave_original.dat'
    if save_fname2 is None:
        fdir,_=os.path.split(dat_fname)
        save_fname2=fdir+'/'+'depth_per_10min.dat'        
        save_fname3=fdir+'/'+'depth_original.dat'
    np.savetxt(save_fname1,S_ifft_wave,fmt='% 1.7f',delimiter='    ')
    np.savetxt(save_fname2,np.mean(dep[:-n].reshape(-1,Fs*time_group),axis=1),fmt='% 1.7f',delimiter='    ')
    np.savetxt(save_fname3,dep,fmt='% 1.7f',delimiter='    ')
    
    if ireturn:
        return S_ifft_wave,dep