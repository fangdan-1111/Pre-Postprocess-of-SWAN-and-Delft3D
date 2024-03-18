import pandas as pd
import numpy as np
import datetime
import shutil
import math

def change_bnd_endpoint(bnd_fname,new_bnd_fname=None):
    # usage:df=change_bnd_endpoint(bnd_fname)
    # usage: 先使用delft3d中的boundary模块规定边界和边界名称，生成一个初始的bnt文件，然后在这里输入文件名，
    # 就会修改文件中的endpoint名称，并在原bnt路径下生成一个后缀为_backup的bnt文件。
    # !!!注意：若边界条件类型不是n个首尾相连的水位边界+1个流量边界，则函数需要再调整
    colspecs = [(0,19),(21,22),(23,24),(25,30),(31,36),(37,42),(43,48),(49,64),(65,77),(78,82)]
    df = pd.read_fwf(bnd_fname,colspecs=colspecs,header=None,delimiter=None)
    df[8][0]=df[0][0]*2
    df[9][0]=df[0][0]+df[0][1]
    df[9][len(df[0])-1]=np.nan
    for i in range(len(df[0])-2):
        j=i+1
        df[9][j]=df[0][j]+df[0][j+1]
        df[8][j]=df[9][i]
    df[9][len(df[0])-2]=df[0].values[-2]*2
    if new_bnd_fname==None:
        new_bnd_fname=bnd_fname[:-4]+'_backup.bnd'
    shutil.copy(bnd_fname,new_bnd_fname)
    lines=[None for _ in range(len(df[0]))]
    for i in range(len(df[0])):
        lines[i]=f'{df[0][i]:21s}{df[1][i]:2s}{df[2][i]:2s}{df[3][i]:5d}{df[4][i]:6d}{df[5][i]:6d}{df[6][i]:6d}  {df[7][i]:13.7e}0 {df[8][i]:13s}{df[9][i]}           \n'
    lines[-1]=lines[-1][0:75]+'          \n'
    with open(bnd_fname, 'w') as f:
        f.writelines(lines)
    return df
    
def change_bnd_index(old_bnd_fname,x_index=None,y_index=None,new_bnd_fname=None):
    # old_bnd_fname=r"E:\d3d_cases\dk_chongming\dept2_caolv35_202212\2D2022_cjk_cmmi5.bnd"
    # new_bnd_fname=r"E:\d3d_cases\chongming_flow\4_bound\chongming_flow_v1.bnd"
    
    import pandas as pd
    import warnings
    import numpy as np
    warnings.filterwarnings("ignore")
    colspecs = [(0,19),(21,22),(23,24),(25,30),(31,36),(37,42),(43,48),(49,64),(65,77),(78,82)]
    df = pd.read_fwf(old_bnd_fname,colspecs=colspecs,header=None,delimiter=None)
    df[5].values[-2]=x_index[-1]
    df[6].values[-2]=1
    if x_index is None:
        dic = read_bnd_nodes(old_bnd_fname)
        old_grid_fname = input("请输入已调好的模型的grd文件'：")
        new_grid_fname = input("请输入本模型的grd文件'：")
        x_coord, y_coord = read_grid_xy(old_grid_fname,dic['x_index'],dic['y_index'])
        x_index, y_index = read_grid_mn(new_grid_fname,x_coord,y_coord)
        x_index=np.array(x_index)+2
        y_index=np.array(y_index)+2
    if new_bnd_fname is None:
        new_bnd_fname=old_bnd_fname[-4:]+'_new.bnd'
    # 根据这里的index写入新的bnd文件

    n_ind=[]
    e_ind=[]
    s_ind=[]
    for i in range(len(x_index)-1):
        bound_name=df[0][i]
        if bound_name[0]=='n' or bound_name[0]=='N':
            n_ind.append(i)
        elif bound_name[0]=='e' or bound_name[0]=='E':
            e_ind.append(i)
        elif bound_name[0]=='s' or bound_name[0]=='S':
            s_ind.append(i)
    # max_m=max(x_index)
    # max_n=max(y_index)
    # min_m=1
    # min_n=1

    for i in range(len(x_index)-1):
        df[3][i]=x_index[i]
        df[4][i]=y_index[i]
    for i in range(len(x_index)-2):
        if i<max(n_ind):
            df[5][i]=df[3][i+1]-1
            df[6][i]=df[4][i+1]
        elif i==max(n_ind):
            df[5][i]=df[3][i+1]-1
            df[6][i]=df[4][i+1]
            df[4][i+1]=df[4][i+1]-1
        elif i<max(e_ind):
            df[5][i]=df[3][max(e_ind)]
            df[6][i]=df[4][i+1]
        elif i==max(e_ind):
            df[4][i+1]=1
            df[5][i]=df[3][i+1]
            df[6][i]=df[4][i+1]+1
            df[3][i+1]=df[3][i+1]-1
        elif i<=max(s_ind):
            df[4][i+1]=1
            df[5][i]=df[3][i+1]+1
            df[6][i]=df[4][i+1]
    
    lines=[None for _ in range(len(df[0]))]
    for i in range(len(df[0])):
        lines[i]=f'{df[0][i]:21s}{df[1][i]:2s}{df[2][i]:2s}{df[3][i]:5d}{df[4][i]:6d}{df[5][i]:6d}{df[6][i]:6d}  {df[7][i]:13.7e}0 {df[8][i]:13s}{df[9][i]}           \n'
    lines[-1]=lines[-1][0:75]+'          \n'
    with open(new_bnd_fname, 'w') as f:
        f.writelines(lines)
    
def read_bnd_nodes(bnd_fname):
    # usage:dic=read_bnd_nodes(new_bnd_fname)
    # 此函数用于读出bnd文件中的边界endpoint名称和结构网格中边界节点的xy—index
    # 注意：不读最后一行是因为最后一行是流量时间序列边界，若边界条件类型不是n个水位边界+1个流量边界则函数需要再调整
    colspecs = [(0,19),(21,22),(23,24),(25,30),(31,36),(37,42),(43,48),(49,64),(65,77),(78,82)]
    df = pd.read_fwf(bnd_fname,colspecs=colspecs,header=None,delimiter=None)
    eps=df[8].values
    eps[-1]=df[9][len(df[0])-2]
    x1=df[3].values[:-1:]
    y1=df[4].values[:-1:]
    x2=df[5].values[:-1:]
    y2=df[6].values[:-1:]
    x=df[3].values
    x[-1]=df[5].values[-2]
    y=df[4].values
    y[-1]=df[6].values[-2]
    dic={'end_points':eps,'x_index1':x1,'y_index1':y1,'x_index2':x2,'y_index2':y2,'x_index':x,'y_index':y}
    return dic

def read_grid_xy(grid_fname,x_index=[0,1],y_index=[0,1],if_allgrd=False):
    '''
        此函数用于读出d3d结构网格.grd文件中，索引为x_index,y_index的节点的xy坐标值
    '''
    with open(grid_fname, 'r') as f:
        lines = [line.rstrip('\n').split() for line in f]
    nx=int(lines[6][0])
    ny=int(lines[6][1])
    x_index=[x-2 if (x-1)>=nx else x-1 for x in x_index]# 问题：grd文件显示nx=2516,ny=252
    y_index=[y-2 if (y-1)>=ny else y-1 for y in y_index]# 但是bnt文件中指示边界格点的xy索引都超出了1格，怎么回事?
    del lines[0:8]#ETA从第m行，开始这里就写0:m-1
    for line in lines:
        if len(line)>5:
            del line[0:2]
    len_each_eta=math.ceil(nx/5)
    all_grd=[None,None]
    for x_or_y in range(2):
        all_grd[x_or_y]=[None for _ in range(ny)]
        for k in range(ny):
            temp=[]
            for i in range(len_each_eta):
                ind=(x_or_y*ny+k)*len_each_eta+i
                for j in range(len(lines[ind])):
                    temp.append(lines[ind][j])
            all_grd[x_or_y][k]=[float(l) for l in temp]
    all_grd=np.array(all_grd)
    n=len(x_index)
    #第1个坐标：想要x坐标还是y坐标
    #第2个坐标：点的y索引
    #第3个坐标：点的x索引
    x=all_grd[np.zeros(n,dtype=int),y_index,x_index]
    y=all_grd[np.ones(n,dtype=int),y_index,x_index]
    if if_allgrd:
        return all_grd
    else:
        return x,y

def read_grid_mn(grid_fname,x_coord=[0,1],y_coord=[0,1]):
    '''
        此函数用于读出d3d结构网格.grd文件中，最接近坐标x_coord，y_coord的节点的m，n索引
    '''
    import numpy as np
    all_grd=read_grid_xy(grid_fname,if_allgrd=True)
    def find_nearest_point(A, all_grd):
        min_distance = float('inf')

        for i in range(all_grd.shape[1]):
            for j in range(all_grd.shape[2]):
                distance = np.sqrt((A[0] - all_grd[0, i, j])**2 + (A[1] - all_grd[1, i, j])**2)
                if distance < min_distance:
                    min_distance = distance
                    y = i
                    x = j
        return x, y
    x_index=[]
    y_index=[]
    for k in range(len(x_coord)):
        A=[x_coord[k],y_coord[k]]
        tem1,tem2=find_nearest_point(A, all_grd)
        x_index.append(tem1)
        y_index.append(tem2)
    return x_index, y_index

def write_bca(dic,bca_fname):
    # 此函数用于创建一个d3d需要的开边界文件（边界类型：天文潮调和常数）
    # 需要输入一个dictionary类型变量，包含'end_points'（values=边界endpoints名称）、
    # 'con_name'（values=考虑的分潮名）、'M2'等'con_name'中包含的分潮名
    # （values=2*len(end_points)的list，列表的第一列是振幅，第二列是相位）
    lines=[]
    for i in range(len(dic['x_index'])):
        lines.append(f"{dic['end_points'][i]}\n")
        for cname in dic['con_name']:
            lines.append(f"{cname} {dic[cname][0][i]:13.7e}  {dic[cname][1][i]:13.7e}\n")
    with open(bca_fname, 'w') as f:
        f.writelines(lines)
        
def write_bct(start_time,end_time,bct_fname,new_bct_fname):
    import sys
    sys.path.append(r"D:\Anaconda\myPythonProj")
    import my_utils

    # 此函数用于创建一个d3d需要的开边界文件（边界类型：时间序列开边界，此处使用的是大通流量边界）
    # 先使用delft3d中的boundary模块生成一个初始的bct文件，然后在这里输入文件名，
    # 就会修改文件中的endpoint名称，并在原bct路径下生成一个后缀为_backup的bct文件。
    # 输入：eg.    start_time=(2021,7,16)  end_time=(2021,7,31)    bct_fname=r"E:\d3d_cases\cjk\cjk_dm1.bct"
    
    # default settings for Datong discharge file
    data_fname=rf'E:\d3d_cases\cjk\dk\{start_time[0]:4d}年大通流量.xlsx'
    sheet_name='Sheet1'
    time_col_name='日期'
    time_input_format="%Y/%m/%d"
    time_output_format="%d %m %Y %H %M %S"
    time_convert=-8# convert from Beijing Time to GMT
    var_col_name='流量（m³/s）'

    st=datetime.datetime(*start_time)
    et=datetime.datetime(*end_time)
    df=pd.read_excel(data_fname,sheet_name=sheet_name)
    tm_str=df[time_col_name].values
    tm_all = [datetime.datetime.now for _ in range(len(tm_str))]
    for i,time in enumerate(tm_str):
        if time.dtype=='datetime64[ns]':
            tm_all[i]=datetime.datetime.fromtimestamp(time.astype(datetime.datetime)/1000000000)+datetime.timedelta(hours=time_convert)
        else:
            tm_all[i]=datetime.datetime.strptime(time, time_input_format)+datetime.timedelta(hours=time_convert)
    ind=my_utils.in_range(tm_all,st,et)
    target_values=np.array(df[var_col_name].values[ind],dtype='float')
    target_times_dt=np.array(tm_all)[ind]

    if new_bct_fname is None:
        new_bct_fname=bct_fname[:-4]+'_new.bct'
    shutil.copy(bct_fname,new_bct_fname)
    with open(bct_fname, 'r') as f:
        lines=f.readlines()
    l10='records-in-table     '+f'{len(target_values)}'+'\n'
    lines[10]=l10
    new_lines=lines[0:11]
    ref_tm=datetime.datetime.strptime(lines[4][-9:-1],'%Y%m%d')
    for i in range(len(target_values)):
        str1=f'{(target_times_dt[i]-ref_tm).total_seconds()/60:13.7e}'
        tm_str=str1[:-1]+'0'+str1[-1]
        dc_str=f'{target_values[i]:=13.7e}'[:-1]+'0'+f'{target_values[i]:13.7e}'[-1]
        end_str='9.9999900e+002'
        new_lines.append(' '+tm_str+'  '+dc_str+'  '+end_str+'\n')
    with open(new_bct_fname, 'w') as f:
        f.writelines(new_lines)
        
def write_rgh(dep_fname,rgh_fname,mm,nn):
    '''
    mm,nn是x、y方向的网格数，dep文件中会有nn+1行,mm+1列
    '''
    with open(dep_fname, 'r') as f:
        lines = [line.rstrip('\n').split() for line in f]
    all_dep=np.array([obj for sublist in lines for obj in sublist],dtype='float32')
    all_dep=all_dep.reshape(nn+1,mm+1)
    n=all_dep**(1/6)/(18*np.log10(12*all_dep/0.005))#糙率公式1（邓珂小论文公式7）
    
    # test：
    n[n>0.022]=0.022
    n[n<0.01]=0.01
    n[np.isnan(n)]=0.022
    n_max=n.max()
    n_min=n[n>0.01].min()
    n=(n-n_min)/(n_max-n_min)*(0.022-0.01)+0.01
    n[n<0.01]=0.01
    
    # 原：
    # n[n>0.022]=0.022
    # n[n<0.01]=0.01
    # n[np.isnan(n)]=0.01
    
    n=np.vstack((n,n))

    # 打开文件以写入模式
    # 此段代码用于将数组写入文本，每行最多12个数字，数组的下一行也会在文本文件中换行
    with open(rgh_fname, "w") as f:
        for row in n:
            count = 0
            for elem in row:
                f.write('{: 1.7E}  '.format(elem)) # 在这里改数字转字符串的格式
                count += 1
                if count == 12: # 在这里改每行最多几个数字
                    f.write('\n')
                    count = 0
            if count > 0:
                f.write('\n')
                
                
def change_obs_index(old_obs_fname,x_index=None,y_index=None,new_obs_fname=None,old_grid_fname=None,new_grid_fname=None):
    # old_obs_fname=r"E:\d3d_cases\dk_chongming\dept2_caolv35_202212\2D2022_cjk06_cmmi.obs"
    # new_obs_fname=r"E:\d3d_cases\chongming_flow\6_outputs\chongming_flow_v1.obs"

    import pandas as pd
    import warnings
    import numpy as np
    warnings.filterwarnings("ignore")
    colspecs = [(0,20),(22,27),(29,34)]
    df = pd.read_fwf(old_obs_fname,colspecs=colspecs,header=None,delimiter=None)

    if x_index is None:
        old_x_ind=[]
        old_y_ind=[]
        for i in range(len(df[0])):
            old_x_ind.append(df[1][i])
            old_y_ind.append(df[2][i])
        if old_grid_fname is None:
            old_grid_fname = input("请输入已调好的模型的grd文件'：")
            new_grid_fname = input("请输入本模型的grd文件'：")
        # old_grid_fname=r"E:\d3d_cases\dk_chongming\dept2_caolv35_202212\allGRIDcm0.5_cmmi.grd"
        # new_grid_fname=r"E:\d3d_cases\chongming_flow\3_grid\chongming_flow_v1.grd"

        x_coord, y_coord = read_grid_xy(old_grid_fname,old_x_ind,old_y_ind)
        x_index, y_index = read_grid_mn(new_grid_fname,x_coord,y_coord)
        x_index=np.array(x_index)+1
        y_index=np.array(y_index)+1
    if new_obs_fname is None:
        new_obs_fname=old_obs_fname[-4:]+'_new.obs'

    for i in range(len(df[0])):
        df[1][i]=x_index[i]
        df[2][i]=y_index[i]
    lines=[None for _ in range(len(df[0]))]
    for i in range(len(df[0])):
        lines[i]=f'{df[0][i]:<22s}{df[1][i]:>4d}  {df[2][i]:>4d}        \n'
    with open(new_obs_fname, 'w') as f:
        f.writelines(lines)
    return x_index, y_index