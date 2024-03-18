'''
    此脚本用于SWAN文件的前后处理
'''
# preprocess
def convert_grd(grid_fname,dep_fname,in_dir='',out_dir=''):
    '''
        此脚本用于将delft3D创建的曲线网格转换为swan需要的网格和水深格式，转换出的网格文件默认保存在d3d网格的相同文件夹下
        fac: SWAN读取网格坐标的精度只有4位小数点，当使用经纬度坐标且网格很密的时候，需要在swn文件读取网格的时候乘以fac以得到较高精度的网格【已测试，没法增加精度】
        eg.:
            grid_fname=r"E:\d3d_cases\test_swan\d3d_grid\curv.grd"# d3d网格文件
            convert_grd(grid_fname)
    '''
    import numpy as np
    import D3Dpre
    import my_utils
    
    new_grd_name=rf"{out_dir}\{grid_fname[0:-3]}swan_grd"#swan网格文件
    new_dep_name = rf"{out_dir}\{dep_fname[0:-3]}swan_dep"
        
    all_grd=D3Dpre.read_grid_xy(rf'{in_dir}\{grid_fname}',if_allgrd=True)
    all_grd[abs(all_grd)>180]=-999.0
    # tmp=all_grd[0,:,:]
    # tmp1=np.hstack((np.vstack((tmp,-999.*np.ones(tmp.shape[1]))),-999.*np.ones(tmp.shape[0]+1).reshape(-1,1)))
    # tmp=all_grd[1,:,:]
    # tmp2=np.hstack((np.vstack((tmp,-999.*np.ones(tmp.shape[1]))),-999.*np.ones(tmp.shape[0]+1).reshape(-1,1)))
    # test = np.concatenate((tmp1,tmp2),axis=0)
    test = np.concatenate((all_grd[0,:,:],all_grd[1,:,:]),axis=0)
    my_utils.write_formated_file(test,new_grd_name)
    # np.savetxt(new_grd_name,test,fmt='% 1.7E',delimiter='    ')
    
    nxny=[all_grd.shape[2],all_grd.shape[1]]
    
    # 以下用于生成网格数（非格点数）对应的swan_dep文件
    with open(rf'{in_dir}\{dep_fname}', 'r') as f:
        lines = [line.rstrip('\n').split() for line in f]
    all_dep=np.array([obj for sublist in lines for obj in sublist],dtype='float32').reshape([nxny[1]+1,nxny[0]+1])[0:-1,0:-1]
    my_utils.write_formated_file(all_dep,new_dep_name)
    # np.savetxt(new_dep_name,all_dep,fmt='% 1.7E',delimiter='    ')
    

def convert_spc(ww3_spc_fname,swan_spc_fname=' ',n_freq_line = 4,n_dir_line = 4,n_spc_per_line = 7,factor = 1.4E-6):
    '''
        此脚本用于将集群中读取的ww3边界谱文件转换为swan谱文件（参考用户手册P47 BOUND SPEC和P137 Appendix D）格式，转换出的网格文件默认保存在ww3文件的相同文件夹下
        # 默认参数的含义：
            n_freq_line = 4 # ww3文件开头的频率坐标有几行
            n_dir_line = 4 # ww3文件开头的方向坐标有几行
            n_spc_per_line = 7 # ww3文件中频谱每行有几个数
            factor = 1.4E-6 # 参考：SWAN用户手册P143
        eg. :
            ww3_spc_fname=r"E:\d3d_cases\test_swan\attribute_files\ww3.21071600.spc""# ww3输出的谱文件
            convert_spc(ww3_spc_fname)
    '''
    import numpy as np
    import D3Dpre
    import shutil 
    import math
    import textwrap
    
    if swan_spc_fname==' ':
        swan_spc_fname = ww3_spc_fname[:-16]+'swan'+ww3_spc_fname[-13:]
    nhed = 1+n_freq_line+n_dir_line
    new_lines='SWAN 1\n$ Data converted by Dan Fang\n$ Project:null   ;   run number:null\nTIME\n  1\n'

    with open(ww3_spc_fname, 'r') as f:
        lines = [line.rstrip('\n').split() for line in f]
    ntable = int((len(lines)-nhed)/(math.ceil(int(lines[0][3])*int(lines[0][4])/n_spc_per_line+1)*int(lines[0][5])+1))
    nline_in_table = math.ceil(int(lines[0][3])*int(lines[0][4])/n_spc_per_line+1)*int(lines[0][5])+1
    nline_in_loc_table = math.ceil(int(lines[0][3])*int(lines[0][4])/n_spc_per_line)+1

    ########先用一个单独的循环提取位置LONLAT，并且把所有头信息写入new_lines###########
    lonlat = []
    table = lines[nhed:(nhed+nline_in_table)]
    for j in range(int(lines[0][5])):
        line_index1 = 1+j*nline_in_loc_table
        loc_table = table[line_index1:(line_index1+nline_in_loc_table)]
        lonlat.append("  "+loc_table[0][3]+ " " +loc_table[0][2]+"\n")
    new_lines+="LONLAT\n  "+str(len(lonlat))+"\n"
    for lonlat_i in lonlat:
        new_lines+=lonlat_i
    new_lines+="RFREQ\n"+lines[0][3]+"\n"
    tem = ["  "+str(float(string2))+"\n" for string1 in lines[1:(n_freq_line+1)] for string2 in string1]
    tem_freq = np.array([float(string2) for string1 in lines[1:(n_freq_line+1)] for string2 in string1])
    for tem1 in tem:
        new_lines+=tem1
    new_lines+="CDIR\n"+lines[0][4]+"\n"
    tem = np.array(sorted(np.array([float(string2) for string1 in lines[(n_freq_line+1):nhed] for string2 in string1])*180/np.pi)).astype('str').tolist()
    tem_cir = np.array([float(string2) for string1 in lines[(n_freq_line+1):nhed] for string2 in string1])*180/np.pi
    for tem1 in tem:
        new_lines+="  "+tem1+"\n"
    new_lines+="QUANT\n  1\nVaDens\nm2/Hz/degr\n  -0.9900E+02\n"
    #############################################################################

    #######################再用一个大循环提取和写入频谱spc#########################
    for i in range(ntable): # 每个时间步
        line_index = nhed+i*nline_in_table
        table = lines[line_index:(line_index+nline_in_table)]
        new_lines+=table[0][0]+'.'+table[0][1]+"\n"
        new_lines+="FACTOR\n"
        new_lines+=f"  {factor:E}\n"
        
        for j in range(int(lines[0][5])): # 每个位置
            
            line_index1 = 1+j*nline_in_loc_table
            loc_table = table[line_index1:(line_index1+nline_in_loc_table)]
            spc = (np.array([float(string2) for string1 in loc_table[1:] for string2 in string1])/factor).astype(int)
            res = sorted(zip(tem_cir,spc.reshape(int(lines[0][3]),int(lines[0][4]))))
            spc = np.array([x for _,x in res]).reshape(-1)
            text = ''
            for tem in spc:
                text+=f"{tem:<5d} "
            new_lines+=textwrap.fill(text ,width=72)#, placeholder='  '
            new_lines+="\n"
    with open(swan_spc_fname, 'w') as f:
        f.writelines(new_lines)
    for orient in ['N','S','W','E']:
        shutil.copy(swan_spc_fname,rf"D:\SWAN\swan\swan.21071600{orient}.spc")

def generate_plant_rf(grid_fname=r"E:\d3d_cases\test_swan\d3d_grid\curv.grd",
                        samples_xlsx=r"E:\OneDrive\OneDrive - stu.ouc.edu.cn\document\13_我的数据\20231012 block3测量\20231012 block3测量.xlsx",
                        samples_veg_edge_csv=r'X:\GEOSTAT\藨草区域20230105岸线.csv',
                        veg_edge_fname = r"E:\d3d_cases\test_swan\attribute_files\整个研究区域的岸线.csv",
                        sirpus_block_fname = r"E:\d3d_cases\test_swan\attribute_files\藨草区域.csv",
                        reed_block_fname = r"E:\d3d_cases\test_swan\attribute_files\芦苇区域.csv",
                        sim_block_fname = r"E:\d3d_cases\test_swan\attribute_files\研究区域减去芦苇.csv",
                        density_fac = 100,isamplegrid=False,ispherical_gird=True,
                        random_seed=20240116,iplot=True):
    '''
        此脚本用于生成植被参数（密度、高度、直径）随机场文件，其中暂时只有密度场可以输入到SWAN中。
        # 参数含义：
        grid_fname delft3d曲线网格文件grid_fname，注意必须是曲线网格
        samples_xlsx 植被样本文件samples_xlsx，格式与r"E:\OneDrive\OneDrive - stu.ouc.edu.cn\document\13_我的数据\20231012 block3测量\20231012 block3测量.xlsx",sheet_name="植被参数"相同
        samples_veg_edge_csv 植被样本区域的岸线文件samples_veg_edge_csv（可以与整个研究区域的岸线相同），格式参考'X:/GEOSTAT/藨草区域20230105岸线.csv'
        coord_type坐标的类型，目前只能是'spherical'【未来需要改进代码和文件，使其也可以是'cartesian'】
        sirpus_block_fname,reed_block_fname,sim_block_fname
            由奥维地图导出的csv文件（注意导出的时候只选经纬度），也可以手动创建，格式参考默认文件。
            在奥维中绘制岸线等对象的时候一定要全部沿逆时针方向绘制！
            换了研究区域的话这几个文件都需要重新制作，如果只换网格没换研究区域的话可以沿用这个文件。
            !!!!!!!!!!!!!!!!!!!!!!!这里必须注意(第四步 line378)!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!因为没有芦苇样本，所以假设整个潮滩都长的是藨草，因此合并了藨草和芦苇两个区域作为新的藨草区域!!!
        9 density_fac表示植被密度应该乘以多少，才能将单位换算成平方米，例如：采样的样方为0.1m*0.1m，那么fac应当设为100：
            100=1/0.01；如果采样样方为0.2m*0.2m，那么fac=25。注意：swn文件中fac为1。
        10 isamplegrid表示是否绘制一个规则网格作为示例
        11 ispherical_gird表示网格是地理坐标还是笛卡尔坐标。如果半方差模型拟合不好的话，调整line174的scale_length和lag。
        
        eg. :
            grid_fname=r"E:\d3d_cases\test_swan\d3d_grid\curv.grd"# d3d网格文件
            generate_plant_rf(grid_fname=grid_fname)
    '''
    import numpy as np
    from matplotlib import pyplot as plt
    from matplotlib.pyplot import get_cmap
    from matplotlib import rcParams
    import matplotlib
    import pandas as pd
    import gstools as gs
    import textwrap
    from scipy.optimize import fmin_cobyla 
    import math
    import my_utils
    import D3Dpre
    import re
    
    colors=my_utils.my_colors('gradual_change',n_colors=10,color_set='green2brown')
    rcParams['font.family'] = 'Times New Roman'
    rcParams['font.size'] = 15
    
    # 第一步：读取植被样本
    # 读取文件并生成gs.field类型数据，如果植被样本文件格式不一样，应修改这里的代码
    df = pd.read_excel(samples_xlsx,sheet_name="Sheet1")
    data_array = np.array(df.set_index('编号'))
    lat=data_array[:,1]
    lon=data_array[:,2]
    x=data_array[:,3]
    y=data_array[:,4]
    ds=data_array[:,5]*density_fac/100
    ht=data_array[:,6]
    ht[0]=np.mean(ht[1::])
    dm=data_array[:,7]
    
    if ispherical_gird:
        coord_unit='[°]'
        scale_length=0.0011
        lag=0.0001
        
    else:
        lon=x
        lat=y
        coord_unit='[m]'
        scale_length=110
        lag=12.35
    
    # 第二步：拟合趋势（这里拟合的是一维趋势）
    # 读取岸线
    
    land_bound=pd.read_csv(samples_veg_edge_csv,header=None).values
    
    #计算样本点和海岸线的最短距离
    target_points=np.array([lon,lat]).T
    sampledis_to_csl=np.array([my_utils.calculate_min_dis(target_points[i,:],land_bound) for i in range(target_points.shape[0])])
    def plt_linear_trend(x,y,ylabel,xlabel='Distance to vegetation edge '+coord_unit,colors=[],ax=[],iplot=True):
        from scipy.stats import t
        if len(colors)==0:
            colors=my_utils.my_colors('contrast',n_colors=10)
        # 设置字体和字号
        # rcParams['font.family'] = 'Times New Roman'
        # plt.rcParams['font.family'] = 'Times New Roman'
        # plt.rcParams['font.size'] = 15
        # 使用 polyfit() 函数进行线性拟合
        coefficients = np.polyfit(x, y, 1)
        slope = coefficients[0]
        intercept = coefficients[1]
        # 计算拟合值
        predicted = slope * x + intercept
        # 计算R2
        r_squared = 1 - (np.sum((y - predicted)**2) / np.sum((y - np.mean(y))**2))
        
        # 得到斜率和截距的标准误差
        n = len(x)
        std_error_slope = np.sqrt(np.sum((y - predicted)**2) / (n - 2)) / np.sqrt(np.sum((x - np.mean(x))**2))
        std_error_intercept = std_error_slope * np.sqrt(np.sum(x**2) / n)
        # 计算置信区间（以95%置信水平为例）
        alpha = 0.05  # 置信水平为95%
        t_critical = t.ppf(1 - alpha/2, n - 2)
        slope_lower, slope_upper = slope - t_critical * std_error_slope, slope + t_critical * std_error_slope
        intercept_lower, intercept_upper = intercept - t_critical * std_error_intercept, intercept + t_critical * std_error_intercept
        if iplot:
            print(f"拟合的R平方值：{r_squared}")
            print(f"斜率置信区间：({slope_lower}, {slope_upper})")
            print(f"截距置信区间：({intercept_lower}, {intercept_upper})")
        # # 绘制原始数据和拟合直线以及置信区间
        # if iplot:
            if ax==[]:
                ax = plt.gca()
            ax.scatter(x, y, label='samples',marker='s',color=colors[0])
            ax.plot(x, slope * x + intercept, color=colors[3], label='fitted')
            xx = np.linspace(0.95*min(x), 1.05*max(x), 10)
            ax.fill_between(xx, (slope_lower * xx + intercept_lower), (slope_upper * xx + intercept_upper), color=colors[0], alpha=0.2, label='95% CI')
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.legend()
            ax.text(0.13,0.82,f'R$^2$ :  {r_squared:5f}',transform=ax.transAxes)
            # plt.savefig(ylabel+'.png',dpi=200)
            # plt.show()
        return slope, intercept
    slope = np.zeros(3)
    intercept = np.zeros(3)
    if iplot:
        fig, axes = plt.subplots(1, 3, figsize=[18,5.5],tight_layout=True)
    else:
        axes=[[],[],[]]

    slope[0],intercept[0] = plt_linear_trend(sampledis_to_csl,ds,"density [stems/0.01m$^2$]",ax=axes[0],iplot=iplot)#,colors=colors)
    slope[1],intercept[1] = plt_linear_trend(sampledis_to_csl,ht,"height [cm]",ax=axes[1],iplot=iplot)#,colors=colors)
    slope[2],intercept[2] = plt_linear_trend(sampledis_to_csl,dm,"diameter [mm]",ax=axes[2],iplot=iplot)#,colors=colors)
    if iplot:
        plt.savefig('3par.png',dpi=200)
        plt.show()
    
    # 第三步：拟合半方差模型，以°为距离单位
    
    def trend(x,y,coast_line,slope,intercept):
        target_points=np.array([x,y]).T
        dis=np.array([my_utils.calculate_min_dis(target_points[i,:],coast_line) for i in range(target_points.shape[0])])
        trend1=dis*slope+intercept
        return trend1
    def read_ov_csv(ovcsv_fname):
        with open(ovcsv_fname, 'r') as f:
            lines = [line.rstrip('\n').split() for line in f]
        return np.array(re.split('[,;|\t]',lines[1][0][1:-1]),dtype="float").reshape(-1,2)
    veg_edge=read_ov_csv(veg_edge_fname)
    ds_resi=ds-np.mean(ds)
    ht_resi=ht-trend(lon,lat,veg_edge,slope[1],intercept[1])
    dm_resi=dm-np.mean(dm)

    ################# 循环定义变量##################
    for i,var_name in enumerate(['ds','ht','dm','ds_resi','ht_resi','dm_resi']):
        exec("%s_field=gs.field.Field(dim=2)"%var_name)
        exec("%s_field(pos=(x,y),field=%s)"%(var_name,var_name))

    def test_vario_model(xy,field,max_lag,bin_wids,models,ishow=True,ax=[],colors=my_utils.my_colors('contrast',n_colors=10),linestyles=['-','--'],return_ax=False,return_bin=False,return_emp_var=False):
        
        scores = {}
        paras = {}
        pcovs = {}
        fit_models = {}
        
        bins = np.arange(0, max_lag, bin_wids)
        bin_center, vario = gs.vario_estimate(
            *(xy, field, bins),
            mesh_type="unstructured"
        )
        if ishow:
            if ax==[]:
                ax = plt.gca()
            ax.scatter(bin_center, vario, color=colors[9], label="empirical",marker='s')
            ax.legend(loc="lower right")
            
            
        for i,model in enumerate(models):
            fit_models[model] = models[model](dim=2)
            para, pcov, r2 = fit_models[model].fit_variogram(bin_center, vario, return_r2=True)
            scores[model] = r2
            paras[model] = para
            pcovs[model] = pcov
            if ishow:
                fit_models[model].plot(x_max=max_lag, ax=ax,color=colors[i],linestyle=linestyles[i//(len(models)//2)])
        ranking = sorted(scores.items(), key=lambda item: item[1], reverse=True)
        if ishow:
            print("RANKING by Pseudo-r2 score")
            for i, (model, score) in enumerate(ranking, 1):
                print(f"{i:>6}. {model:>15}: {score:.5}")
        best_model_type = ranking[0][0]
        best_model = fit_models[best_model_type]

        return paras[best_model_type], pcovs[best_model_type], scores, best_model,ax,bin_center,vario

    ######绘制拟合半方差模型的图像
    models = {
        "Gaussian": gs.Gaussian,
        "Exponential": gs.Exponential,
        "Matern": gs.Matern,
        "Stable": gs.Stable,
        "Rational": gs.Rational,
        "Circular": gs.Circular,
        "Spherical": gs.Spherical,
        "SuperSpherical": gs.SuperSpherical,
        #"JBessel": gs.JBessel,
    }
    # variables=['ds','ht','dm','ds_resi','ht_resi','dm_resi']
    variables=['ds_resi','ht_resi','dm_resi']
    units=[' [stems/m$^2$]',' [cm]',' [mm]',' [stems/0.01m$^2$]',' [cm]',' [mm]']
    var_long_name=['density','height','diamter','density','height','diamter']

    # paras={}
    # pcovs={}
    r2s={}
    bst_mds={}
    if iplot:
        rcParams['font.size'] = 7
        fig, axes = plt.subplots(1, 3, figsize=[12, 3],tight_layout=True)

    for i,var in enumerate(variables):
        _,_,r2s[var],bst_mds[var] ,_,_,_= test_vario_model((lon,lat),eval(var+'_field'),scale_length,lag,models,ax=axes[np.mod(i,3)],ishow=iplot)#i//3,
        if iplot:
            axes[np.mod(i,3)].set_title('fitting plant '+var_long_name[i]+' to an isotropic model')
            axes[np.mod(i,3)].set_xlabel('distance '+coord_unit)
            axes[np.mod(i,3)].set_ylabel('semivariance')
            txt='        R$^2$: \n'
            for model,r2 in r2s[var].items():
                txt=txt+f'{model:>15}: {r2:<.5}\n'
            axes[np.mod(i,3)].text(0.2,0,txt,transform=axes[np.mod(i,3)].transAxes, ha='left')
    if iplot:
        plt.savefig('compare_vario_models.png',dpi=200)
    
    # 第四步：读取网格、对应的岸线、藨草和芦苇区域
    # initialize
    if not isamplegrid:
        # nplants_fname=grid_fname[0:-3]+"swan_nplants" # swan植被密度场文件名 idla=4 nhedf=0
        # pla_fname=grid_fname[0:-3]+pla # d3d植被密度场文件

        all_grd=D3Dpre.read_grid_xy(grid_fname,if_allgrd=True)
        nxny=[all_grd.shape[2],all_grd.shape[1]]
        points = all_grd.reshape([2,-1]).T
        all_X = points[:,0]
        all_Y = points[:,1]

        # 读取岸线.csv、藨草区域.csv、芦苇区域.csv：
        veg_edge=read_ov_csv(veg_edge_fname)
        sirpus_block=read_ov_csv(sirpus_block_fname)
        reed_block=read_ov_csv(reed_block_fname)
        sim_block=read_ov_csv(sim_block_fname)
        #注意这里！因为没有芦苇样本，所以假设整个潮滩都长的是藨草，因此合并了藨草和芦苇两个区域作为新的藨草区域:
        new_sirpus_block = my_utils.get_new_block(sirpus_block,reed_block)

        from shapely import geometry

        points = all_grd.reshape([2,-1]).T
        polygon = geometry.Polygon(new_sirpus_block)
        res=[]
        res.append([polygon.covers(geometry.Point(point)) for point in points])
        res = np.array(res).squeeze()
        X = points[res,:][:,0]
        Y = points[res,:][:,1]

    # 第五步：生成随机场并绘图
    # 设置字体字号
    rcParams['font.family'] = 'Times New Roman'
    rcParams['font.size'] = 12
    # plt.rcParams['font.family'] = 'Times New Roman'
    # plt.rcParams['font.size'] = 12

    paras={}
    pcovs={}
    r2s={}
    bst_mds={}
    bins={}
    varios={}
    if isamplegrid:
        # x_grid=np.linspace(min(lon)*0.999985,max(lon)*1.000015,50) # 不绘制样本点
        # y_grid=np.linspace(min(lat)*0.99995,max(lat)*1.00005,50)
        # [X1,Y1]=np.meshgrid(x_grid-min(x_grid),y_grid-min(y_grid))
        # X=X1.reshape(-1)
        # Y=Y1.reshape(-1)
        def calculate_which_hand(target_point,lines):
            # 计算二维空间内点到线段的最短距离，lines的格式是2*
            import math
            def check_list_and_convert(variable):
                if isinstance(variable, list):
                    variable = np.array(variable)
                return variable
            check_list_and_convert(target_point)
            check_list_and_convert(lines)    
            def points_to_segments(points):
                segments = []
                for i in range(points.shape[0] - 1):
                    segment = [*points[i], *points[i + 1]]
                    segments.append(segment)
                return np.array(segments)
            def point_to_segment_distance(x, y, x1, y1, x2, y2):
                def distance(x1, y1, x2, y2):
                    return np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
                def dot_product(v1, v2):
                    return v1[0] * v2[0] + v1[1] * v2[1]
                def vector_length(v):
                    return np.sqrt(v[0]**2 + v[1]**2)
                px = x
                py = y
                x1y1 = np.array([x1, y1])
                x2y2 = np.array([x2, y2])
                if x1 == x2 and y1 == y2:
                    return distance(px, py, x1, y1)
                line_vec = x2y2 - x1y1
                point_vec = np.array([px - x1, py - y1])
                line_length = vector_length(line_vec)
                line_unitvec = line_vec / line_length
                dot = dot_product(point_vec, line_unitvec)
                if dot < 0:
                    return distance(px, py, x1, y1)
                if dot > line_length:
                    return distance(px, py, x2, y2)
                perpendicular = x1y1 + line_unitvec * dot
                which_hand = math.sign(np.cross(np.array([px-x1,py-y1]),np.array([x2-x1,y2-y1])))
                return distance(px, py, perpendicular[0], perpendicular[1]),which_hand
            # 初始化最小距离为一个很大的数
            min_distance = float('inf')
            # 将点集合转换为线段集合
            segments=points_to_segments(lines)
            # 遍历每个线段，计算最短距离
            distances=np.zeros(segments.shape[0])
            which_hands=np.zeros(segments.shape[0])
            for i,segment in enumerate(segments):
                x1, y1, x2, y2 = segment
                distances[i],which_hands[i] = point_to_segment_distance(target_point[0], target_point[1], x1, y1, x2, y2)
            min_distance = min(distances)
            which_hand = which_hands[distances==min_distance]
            return which_hand
        
        x_grid=np.linspace(min(lon)*0.999985,max(lon)*1.000015,50) # 绘制样本点
        y_grid=np.linspace(min(lat)*0.99995,max(lat)*1.00005,50)
        [X1,Y1]=np.meshgrid(x_grid,y_grid)
        X=X1.reshape(-1)
        Y=Y1.reshape(-1)
        signs=np.array([calculate_which_hand(X[i],y[i],veg_edge) for i in range(len(X))])

    for i,var in enumerate(variables):
        # i = 1                  # debugging
        # var = variables[i]     # debugging
        paras[var],pcovs[var],r2s[var],bst_mds[var],_,bins[var],varios[var] = test_vario_model((lon,lat),eval(var+'_field'),scale_length,lag,models,ishow=False)
        model = bst_mds[var]
        if iplot:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=[9,4],tight_layout=True)
            ax1.scatter(bins[var],varios[var], label="empirical",color=colors[0],marker='s')
            model.plot("vario_axis", ax=ax1, x_max=0.0012, label="model",color=colors[1],linestyle='--')
            txt=textwrap.fill(str(model),width=30, placeholder='   ')
            ax1.text(0.36,0.3,txt,transform=ax1.transAxes)

            ax1.legend(loc="lower right")
            ax1.set_title('fitting plant '+var_long_name[i]+' to an isotropic model')
            ax1.set_xlabel("distance "+coord_unit)
            ax1.set_ylabel("semivariance")
        srf = gs.SRF(model, seed=random_seed)

        field = srf((X, Y))
        # srf.plot(ax=ax2)
        # ax2.set_title('isotropic gaussian random field')
        # ax2.set_aspect("equal")

        srf1=srf
        if var=='ds_resi':
            srf1.field=(field+np.mean(ds))*100
            minvar = 0
        elif var=='ht_resi':
            minvar = 20
            maxvar = 85
            if isamplegrid:
                srf1.field=field+slope[i]*Y+intercept[i]
                srf1.field[srf1.field>maxvar]=maxvar    
            else:
                srf1.field=field+trend(X,Y,veg_edge,slope[i],intercept[i])
                srf1.field[srf1.field>maxvar]=maxvar
        elif var=='dm_resi':
            minvar = 0
            srf1.field=field+np.mean(dm)
        srf1.field[srf1.field<minvar]=0

        norm = matplotlib.colors.Normalize(vmin=max(minvar,min(srf1.field)), vmax=max(srf1.field))
        if iplot:
            if isamplegrid:
                srf1.field=signs*srf1.field
                tcf=ax2.tricontourf(X,Y,srf1.field,cmap=my_utils.my_colors('colormap'),norm=norm)
                # sct=ax2.scatter(lon-min(x_grid),lat-min(y_grid),c=eval(var[0:2]+'_field.field'),cmap=my_utils.my_colors('colormap'),norm=norm,edgecolor='gray',label='samples',linewidth=0.2)
                ax2.legend(bbox_to_anchor=(1.05, 0.95),loc=4,borderaxespad=0.5,frameon=False)#facecolor=None,edgecolor=None,
                ax2.set_aspect("equal")
                
            if not isamplegrid:
                tcf = ax2.scatter(X,Y,c=srf1.field,cmap=my_utils.my_colors('colormap'),norm=norm,s=0.5)
                ax2.fill(sim_block[:,0],sim_block[:,1],color=colors[0],alpha=0.1,edgecolor = 'gray')
                ax2.text(0.1,0.9,'SCIRPUS',transform=ax2.transAxes)
                ax2.text(0.4,0.4,'OCEAN',transform=ax2.transAxes)
            cb=plt.colorbar(tcf)
            cb.ax.set_ylabel(var_long_name[i]+units[i])
            ax2.set_title(var_long_name[i]+" in each grid")
            ax2.set_ylabel("latitude "+coord_unit)
            ax2.set_xlabel("longitude "+coord_unit)
            
            plt.savefig(var_long_name[i]+'.png',dpi=200)
        
        #最后一步：写入文件 #idla=4
        if not isamplegrid:
            if i == 0: # 如果需要生成直径和高度的随机场，注释这一行并调整缩进
                var_in_each_grid = np.zeros(len(points))
                var_in_each_grid[res]=srf1.field
                # tem = var_in_each_grid.reshape(nxny[0],nxny[1])
                # tem1 = np.vstack((np.zeros(nxny[1]),tem))
                # var_in_each_grid = np.hstack((tem1,np.expand_dims(np.zeros(nxny[0]+1),axis=1))).reshape(-1)
                lines=[]
                for j in range(len(var_in_each_grid)):
                    if (j+1)%10 == 0:
                        lines.append(f'{var_in_each_grid[j]:<.3f} \n')
                    else:
                        lines.append(f'{var_in_each_grid[j]:<.3f} ')
                with open(f'{var_long_name[i]}_rs{str(random_seed).zfill(8)[-4:]}.swan_file','w') as f:
                    f.writelines(lines)
                    
                # var_in_each_grid1 = np.zeros(len(points))  # 如果需要生成常数场，取消注释以下行
                # facs = [100,1,1]
                # var_in_each_grid1[res]=eval(var[0:2]).mean()*facs[i]            
                # lines=[]
                # for j in range(len(var_in_each_grid1)):
                    # if (j+1)%10 == 0:
                        # lines.append(f'{var_in_each_grid1[j]:<.3f} \n')
                    # else:
                        # lines.append(f'{var_in_each_grid1[j]:<.3f} ')
                # with open(var_long_name[i]+"_constant.swan_file",'w') as f:
                    # f.writelines(lines)

def generate_nplant_file(grid_fname=r"E:\d3d_cases\test_swan\d3d_grid\curv.grd",
                        proj_name='swan',
                        samples_xlsx=r"E:\OneDrive\OneDrive - stu.ouc.edu.cn\document\13_我的数据\20231012 block3测量\20231012 block3测量.xlsx",
                        samples_veg_edge_csv=r'X:\GEOSTAT\藨草区域20230105岸线.csv',
                        veg_edge_fname = r"E:\d3d_cases\test_swan\attribute_files\整个研究区域的岸线.csv",
                        sirpus_block_fnames = r"E:\d3d_cases\test_swan\attribute_files\藨草区域.csv",
                        sim_block_fname = r"E:\d3d_cases\test_swan\attribute_files\研究区域减去芦苇.csv",
                        density_fac = 100,isamplegrid=False,ispherical_gird=True,
                        random_seed=20240116,iplot=True,icons=True):
    '''
        此脚本用于生成植被密度文件，包括随机场和均匀场
        # 参数含义：
        grid_fname delft3d曲线网格文件grid_fname，注意必须是曲线网格
        samples_xlsx 植被样本文件samples_xlsx，格式与r"E:\OneDrive\OneDrive - stu.ouc.edu.cn\document\13_我的数据\20231012 block3测量\20231012 block3测量.xlsx",sheet_name="植被参数"相同
        samples_veg_edge_csv 植被样本区域的岸线文件samples_veg_edge_csv（可以与整个研究区域的岸线相同），格式参考'X:/GEOSTAT/藨草区域20230105岸线.csv'
        veg_edge_fname，sirpus_block_fnames，sim_block_fname
            由奥维地图导出的csv文件（注意导出的时候只选经纬度），也可以手动创建，格式参考默认文件。
            在奥维中绘制岸线等对象的时候一定要全部沿逆时针方向绘制！
            换了研究区域的话这几个文件都需要重新制作，如果只换网格没换研究区域的话可以沿用这些文件。
            sirpus_block_fnames可以是文件名字符串，也可以是包含多个文件名的lisy,表示有多个藨草区域。
        density_fac表示植被密度应该乘以多少，才能将单位换算成平方米，例如：采样的样方为0.1m*0.1m，那么fac应当设为100：
            100=1/0.01；如果采样样方为0.2m*0.2m，那么fac=25。注意，swn文件里的fac为1。
        ispherical_gird表示网格是地理坐标还是笛卡尔坐标。如果半方差模型拟合不好的话，调整line174的scale_length和lag。
        
        eg. :
            grid_fname=r"E:\d3d_cases\test_swan\d3d_grid\curv.grd"# d3d网格文件
            generate_plant_rf(grid_fname=grid_fname)
    '''
    import numpy as np
    from matplotlib import pyplot as plt
    from matplotlib.pyplot import get_cmap
    from matplotlib import rcParams
    import matplotlib
    import pandas as pd
    import gstools as gs
    import textwrap
    from scipy.optimize import fmin_cobyla 
    import math
    import my_utils
    import D3Dpre
    import re
    
    colors=my_utils.my_colors('gradual_change',n_colors=10,color_set='green2brown')
    rcParams['font.family'] = 'Times New Roman'
    rcParams['font.size'] = 15
    
    # 第一步：读取植被样本
    # 读取文件并生成gs.field类型数据，如果植被样本文件格式不一样，应修改这里的代码
    df = pd.read_excel(samples_xlsx,sheet_name="Sheet1")
    data_array = np.array(df.set_index('编号'))
    lat=data_array[:,1]
    lon=data_array[:,2]
    x=data_array[:,3]
    y=data_array[:,4]
    ds=data_array[:,5]*density_fac/100
    ht=data_array[:,6]
    ht[0]=np.mean(ht[1::])
    dm=data_array[:,7]
    
    if ispherical_gird:
        coord_unit='[°]'
        scale_length=0.0011
        lag=0.0001
        
    else:
        lon=x
        lat=y
        coord_unit='[m]'
        scale_length=110
        lag=12.35
        
    def read_ov_csv(ovcsv_fname):
        with open(ovcsv_fname, 'r') as f:
            lines = [line.rstrip('\n').split() for line in f]
        return np.array(re.split('[,;|\t]',lines[1][0][1:-1]),dtype="float").reshape(-1,2)
    veg_edge=read_ov_csv(veg_edge_fname)
    ds_resi=ds-np.mean(ds)


    ################# 循环定义变量##################
    for i,var_name in enumerate(['ds','ds_resi']):
        exec("%s_field=gs.field.Field(dim=2)"%var_name)
        exec("%s_field(pos=(x,y),field=%s)"%(var_name,var_name))

    def test_vario_model(xy,field,max_lag,bin_wids,models,ishow=True,ax=[],colors=my_utils.my_colors('contrast',n_colors=10),linestyles=['-','--'],return_ax=False,return_bin=False,return_emp_var=False):
        
        scores = {}
        paras = {}
        pcovs = {}
        fit_models = {}
        
        bins = np.arange(0, max_lag, bin_wids)
        bin_center, vario = gs.vario_estimate(
            *(xy, field, bins),
            mesh_type="unstructured"
        )
        if ishow:
            if ax==[]:
                ax = plt.gca()
            ax.scatter(bin_center, vario, color=colors[9], label="empirical",marker='s')
            ax.legend(loc="lower right")
            
            
        for i,model in enumerate(models):
            fit_models[model] = models[model](dim=2)
            para, pcov, r2 = fit_models[model].fit_variogram(bin_center, vario, return_r2=True)
            scores[model] = r2
            paras[model] = para
            pcovs[model] = pcov
            if ishow:
                fit_models[model].plot(x_max=max_lag, ax=ax,color=colors[i],linestyle=linestyles[i//(len(models)//2)])
        ranking = sorted(scores.items(), key=lambda item: item[1], reverse=True)
        if ishow:
            print("RANKING by Pseudo-r2 score")
            for i, (model, score) in enumerate(ranking, 1):
                print(f"{i:>6}. {model:>15}: {score:.5}")
        best_model_type = ranking[0][0]
        best_model = fit_models[best_model_type]

        return paras[best_model_type], pcovs[best_model_type], scores, best_model,ax,bin_center,vario

    ######绘制拟合半方差模型的图像
    models = {
        "Gaussian": gs.Gaussian,
        "Exponential": gs.Exponential,
        "Matern": gs.Matern,
        "Stable": gs.Stable,
        "Rational": gs.Rational,
        "Circular": gs.Circular,
        "Spherical": gs.Spherical,
        "SuperSpherical": gs.SuperSpherical,
        #"JBessel": gs.JBessel,
    }

    variables=['ds_resi']
    units=[' [stems/m$^2$]',' [stems/0.01m$^2$]']
    var_long_name=['density','density']

    # paras={}
    # pcovs={}
    r2s={}
    bst_mds={}
    if iplot:
        rcParams['font.size'] = 7
        fig, axes = plt.subplots(1, 1, figsize=[4, 3],tight_layout=True)

    for i,var in enumerate(variables):
        _,_,r2s[var],bst_mds[var] ,_,_,_= test_vario_model((lon,lat),eval(var+'_field'),scale_length,lag,models,ax=axes,ishow=iplot)#i//3,
        if iplot:
            axes.set_title('fitting plant '+var_long_name[i]+' to an isotropic model')
            axes.set_xlabel('distance '+coord_unit)
            axes.set_ylabel('semivariance')
            txt='        R$^2$: \n'
            for model,r2 in r2s[var].items():
                txt=txt+f'{model:>15}: {r2:<.5}\n'
            axes.text(0.2,0,txt,transform=axes.transAxes, ha='left')
    if iplot:
        plt.savefig('compare_vario_models.png',dpi=200)
    
    # 第四步：读取网格、对应的岸线、藨草和芦苇区域
    all_grd=D3Dpre.read_grid_xy(grid_fname,if_allgrd=True)
    nxny=[all_grd.shape[2],all_grd.shape[1]]
    points = all_grd.reshape([2,-1]).T
    all_X = points[:,0]
    all_Y = points[:,1]

    # 读取岸线.csv、藨草区域.csv、芦苇区域.csv、模拟区域.csv：
    veg_edge=read_ov_csv(veg_edge_fname)
    sim_block=read_ov_csv(sim_block_fname)
    if not type(sirpus_block_fnames)==str:
        sirpus_blocks=[read_ov_csv(fnm) for fnm in sirpus_block_fnames]
    else:
        sirpus_blocks=[read_ov_csv(sirpus_block_fnames)]

    from shapely import geometry
    points = all_grd.reshape([2,-1]).T

    res1=[]
    for sirpus_block in sirpus_blocks:
        polygon = geometry.Polygon(sirpus_block)
        res1.append([polygon.covers(geometry.Point(point)) for point in points])
        
    res1 = np.array(res1).squeeze()

    res=np.zeros(res1.shape[1],dtype='bool')
    for i in range(res1.shape[0]):
        res=np.ma.mask_or(res1[i],res)

    X = points[res,:][:,0]
    Y = points[res,:][:,1]

    # 第五步：生成随机场并绘图
    # 设置字体字号
    rcParams['font.family'] = 'Times New Roman'
    rcParams['font.size'] = 12
    # plt.rcParams['font.family'] = 'Times New Roman'
    # plt.rcParams['font.size'] = 12

    paras={}
    pcovs={}
    r2s={}
    bst_mds={}
    bins={}
    varios={}

    for i,var in enumerate(variables):
        # i = 1                  # debugging
        # var = variables[i]     # debugging
        paras[var],pcovs[var],r2s[var],bst_mds[var],_,bins[var],varios[var] = test_vario_model((lon,lat),eval(var+'_field'),scale_length,lag,models,ishow=False)
        model = bst_mds[var]
        if iplot:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=[9,4],tight_layout=True)
            ax1.scatter(bins[var],varios[var], label="empirical",color=colors[0],marker='s')
            model.plot("vario_axis", ax=ax1, x_max=0.0012, label="model",color=colors[1],linestyle='--')
            txt=textwrap.fill(str(model),width=30, placeholder='   ')
            ax1.text(0.36,0.3,txt,transform=ax1.transAxes)

            ax1.legend(loc="lower right")
            ax1.set_title('fitting plant '+var_long_name[i]+' to an isotropic model')
            ax1.set_xlabel("distance "+coord_unit)
            ax1.set_ylabel("semivariance")
        srf = gs.SRF(model, seed=random_seed)

        field = srf((X, Y))


        srf1=srf
        if var=='ds_resi':
            srf1.field=(field+np.mean(ds))*100
            minvar = 0

        norm = matplotlib.colors.Normalize(vmin=max(minvar,min(srf1.field)), vmax=max(srf1.field))
        if iplot:
            tcf = ax2.scatter(X,Y,c=srf1.field,cmap=my_utils.my_colors('colormap'),norm=norm,s=0.5)
            ax2.fill(sim_block[:,0],sim_block[:,1],color=colors[0],alpha=0.1,edgecolor = 'gray')
            ax2.text(0.1,0.9,'SCIRPUS',transform=ax2.transAxes)
            ax2.text(0.4,0.4,'OCEAN',transform=ax2.transAxes)
            cb=plt.colorbar(tcf)
            cb.ax.set_ylabel(var_long_name[i]+units[i])
            ax2.set_title(var_long_name[i]+" in each grid")
            ax2.set_ylabel("latitude "+coord_unit)
            ax2.set_xlabel("longitude "+coord_unit)
            
            plt.savefig(proj_name+'.png',dpi=200)
        
        #最后一步：写入文件 #idla=4
        var_in_each_grid = np.zeros(len(points))
        var_in_each_grid[res]=srf1.field

        lines=[]
        for j in range(len(var_in_each_grid)):
            if (j+1)%10 == 0:
                lines.append(f'{var_in_each_grid[j]:<.3f} \n')
            else:
                lines.append(f'{var_in_each_grid[j]:<.3f} ')
        with open(f'{proj_name}_rs{str(random_seed).zfill(8)[-4:]}.swan_pla','w') as f:
            f.writelines(lines)
            
        if icons:  # 生成常数场
            var_in_each_grid1 = np.zeros(len(points))
            facs = [100,1,1]
            var_in_each_grid1[res]=eval(var[0:2]).mean()*facs[i]            
            lines=[]
            for j in range(len(var_in_each_grid1)):
                if (j+1)%10 == 0:
                    lines.append(f'{var_in_each_grid1[j]:<.3f} \n')
                else:
                    lines.append(f'{var_in_each_grid1[j]:<.3f} ')
            with open(proj_name+".constant_swan_pla",'w') as f:
                f.writelines(lines)




# postprocess
def convert_mat_to_nc(mat_fname, out_nc_fname=' '):
    '''
    将SWAN生成的mat文件转为一个nc文件（包括所有时间步），方便读取和绘图
    使用：
        mat_fname = r"D:\SWAN\swan\SWANOUT1.mat"
        out_nc_fname = 'SWANOUT1_converted.nc'
        convert_mat_to_nc(mat_fname, out_nc_fname)
    注意：
        读取生成的nc文件最好用nc.Dataset(out_nc_fname)['变量名如Hsig'][第n时间步,:]读取，这样读出的数据还是np掩码数组
        而xr.open_dataset无法识别掩码数组
        可使用ds.variables.keys()查看所有变量名
    '''
    import re
    import numpy as np
    import datetime
    import netCDF4 as nc
    import spicy
    if out_nc_fname==' ':
        out_nc_fname = mat_fname[:-3]+'nc'
    matfile = spicy.io.loadmat(mat_fname)
    keys = np.array(list(matfile.keys()))
    tem = np.array([sub for sub in [re.findall('__(.*?)__',txt) for txt in keys] if len(sub) > 0])
    comment_keys = np.array([tem2 for tem1 in tem for tem2 in tem1])

    coords_keys = np.array([sub for sub in keys if len(sub) == 2] )

    tem = np.array([sub for sub in [re.findall('_([0-9]{8}_[0-9]{6})',txt) for txt in keys] if len(sub) > 0] )
    times = np.array(sorted([sub for sub in list(set([tem2 for tem1 in tem for tem2 in tem1]))]))
    
    datetimes = sorted([datetime.datetime.strptime(time,'%Y%m%d_%H%M%S') for time in times])
    # datetimes = [datetime.datetime.strptime(time,'%Y%m%d_%H%M%S') for time in times]
    
    tem = [sub for sub in [re.findall('(.+[a-z]{0,2}[0-9]{0,2}_[0-9]{8})',txt) for txt in keys] if len(sub) > 0]
    nonsta_keys = list(set(np.array([tem2[0:-9] for tem1 in tem for tem2 in tem1])))

    ds = nc.Dataset(out_nc_fname,'w',format = "NETCDF4")
    coords_var_name = ['longitudes','latitudes']

    ds.createDimension('lonlat',matfile[coords_keys[0]].shape[0]*matfile[coords_keys[0]].shape[1])
    ds.createDimension('time',len(datetimes))
    ds.createVariable('longitudes',np.float32,('lonlat',))
    ds.createVariable('latitudes',np.float32,('lonlat',))
    ds.createVariable('times','str',('time',))

    for key in nonsta_keys:
        ds.createVariable(key, np.float32, ('time','lonlat'))
    for key in comment_keys:
        ds.setncattr(key,matfile[f'__{key}__'])
    ds.setncattr('nxny',str(matfile[coords_keys[0]].shape[0])+' '+str(matfile[coords_keys[0]].shape[1]))

    ds['times'][:] = np.array([str(np.datetime64(dt)) for dt in datetimes])

    for i,var in enumerate(coords_keys):
        if i == 0:
            mask = np.logical_or(np.isnan(matfile[var]),np.abs(matfile[var])>180)
        ds[coords_var_name[i]][:] = np.ma.array(matfile[var],mask=mask).reshape(-1)

    for var in nonsta_keys:
        for i,time in enumerate(times):
            ori_var_name = f'{var}_{time}'
            ds[var][i,:] = np.ma.array(matfile[ori_var_name],mask=mask).reshape(1,-1)
    ds.close()
    
def swan_plot(nc_fname,timesteps=0,var_names = 'Hsig',units = 'm',ax=[],fig=[],var_lims=None,suffix=''):
    '''
        此脚本用于绘制SWANprocess.convert_mat_to_nc()转换出的nc文件
        参数含义：
            nc_fname：SWANprocess.convert_mat_to_nc()转换出的nc文件名
            timesteps：可以是int，表示绘制第几个时间步的数据；
                      也可以是str(int)如'2'，表示以2个时间步的间隔绘制动图
            var_names：绘制哪个变量。可以用ds.variables.keys()查看nc文件的所有变量名。
                可以是str也可是list，若为list则绘制多个变量
            units：变量的单位
                可以是str也可是list，若为list则绘制多个变量
            ax：图像绘制到哪里，若未指定则新建一空白画布fig和轴ax【若绘制多个变量，这个参数将被忽略】
            var_lims：绘图时变量范围的上下限，绘制动图的时候需要用到，保证每个时间步的colorbar是统一的。
                可以是None，可以是元素为float的长度为2的list（如[1,1.5]），也可以是多个上述list组成的list，对应地绘制多个变量
    '''
    import matplotlib.pyplot as plt
    import numpy as np
    import cartopy.crs as ccrs   
    import cartopy.feature as cfeature
    from matplotlib.pyplot import get_cmap
    import netCDF4 as nc
    import matplotlib.ticker as ticker
    import matplotlib
    import os
    import datetime
    import my_utils
    import math
    from matplotlib import rcParams
    
    rcParams['font.size'] = 7.5

    # initialize
    ds = nc.Dataset(nc_fname,'r')
    nxny = np.array(str.split(ds.nxny),dtype='int')
    times = ds['times'][:]

    def plot1(timestep,var_name,unit,var_lim,formater = 'svg',ax=ax,fig=fig,isave=True,suffix='',out_dir=''):
        variables = ds[var_name][:,:]
        if ax==[]:
            fig,ax = plt.subplots()
        lon = ds['longitudes'][:].reshape(nxny[0],nxny[1])
        lat = ds['latitudes'][:].reshape(nxny[0],nxny[1])
        variable = variables[timestep,:].reshape(nxny[0],nxny[1])
        if not var_lims==None:
            norm = matplotlib.colors.Normalize(vmin=var_lim[0],vmax=var_lim[1])
            levels=list(np.linspace(var_lim[0],var_lim[1],num=10,dtype='float32',endpoint=False))
        else:
            levels=var_lim
            norm = None
        tcf = ax.contourf(lon, lat, variable, 10,
                         cmap=my_utils.my_colors('colormap'),levels=levels,extend='both',norm=norm)#
        # tcf = ax.scatter(lon, lat, c=variable,
                         # cmap=my_utils.my_colors('colormap'),norm=norm)
        
        ax.set_xlim([lon.min(),lon.max()])
        ax.set_ylim([lat.min(),lat.max()])
        ax.set_aspect("equal")
        
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad="5%")
        plt.colorbar(tcf, shrink=.98,cax=cax)
        

        ax.set_title(f"{var_name} [{unit}]\ntime:{times[timestep]}")
        if isave:
            plt.savefig(f"{out_dir}{var_name}{''.join(filter(str.isdigit, times[timestep]))[:-10]}output{suffix}.{formater}",format = formater,dpi=200)
    def plot2(timestep,formater = 'svg',suffix='',out_dir=''):#同时绘制多个变量
        nn = len(var_names)
        m = math.floor(nn**0.5)
        n = math.ceil(nn/m)
        fig,axes =plt.subplots(m,n,tight_layout=True)
        for i in range(len(var_names)):
            if m==1:
                ax=axes[i]
            else:
                ax=axes[i//n,np.mod(i,n)]
            var_lim=var_lims if var_lims==None else var_lims[i]
            plot1(timestep,var_names[i],units[i],var_lim,formater = formater,ax=ax,fig=fig,isave=False)
        plt.savefig(f"{out_dir}{len(var_names)}vars{''.join(filter(str.isdigit, times[timestep]))[:-10]}output{suffix}.{formater}",format = formater,dpi=200)
    if type(timesteps) == int:
        if type(var_names)==str:
            plot1(timesteps,var_names,units,var_lims,ax=ax,fig=fig)
        else:
            plot2(timesteps)
    else:
        if type(var_names)==str:
            var_name=var_names
        else:
            var_name=str(len(var_names))+'vars'
        timesteps=np.arange(len(times))[::int(timesteps)]
        import glob
        import matplotlib.animation as anim
        import imageio

        def convert_img2gif(var_name='1', imgs_dir=''):
            files = glob.glob(f'{imgs_dir}*_temp.png')
            files.sort()
            frames = []
            for file in files:
                frames.append(imageio.imread(file))
            imageio.mimsave(f"{var_name}.gif", frames, 'GIF', duration = 500,loop=0)
            for file in files:
                os.remove(file)

        for i,timestep in enumerate(timesteps):
            if type(var_names)==str:
                plot1(timestep,var_names,units,var_lims,formater='png',suffix='_temp',out_dir='tempimage/')
            else:
                plot2(timestep,formater='png',suffix='_temp',out_dir='tempimage/')
        convert_img2gif(var_name=var_name,imgs_dir='tempimage/')

    ds.close()
    
class Spc:
    '''
    ww3_fname：要读取的ww3 output文件
    得到的Spc类变量具有属性：
        fname
        freqs
        dirs
        locs
        times
        spcs：spcs是4维ndarray，四个维度分别表示时间，位置，频率和方向
    得到的Spc类变量具有函数：
        plot_spc(self,nstep,nloc,ax=None,cmap=None,arrowdir=45)：
            参数 nstep--第几个时间步,nloc--第几个位置,
            ax=None,cmap=None,arrowdir=45--频率（半径）指示轴的方向，角度
            ymax=None--频率（半径）的上限，默认参数为自适应上限
    example：
        ww3_fname=r"E:\OneDrive\OneDrive - stu.ouc.edu.cn\document\12_软件学习\SWAN\一些帮助理解的delft3d自动产生的swan文件\ww3.21010100.spc"
        tem=Spc(ww3_fname)
        tem.dirs
        tem.plot_spc(0,0,arrowdir=22.5)
    '''
    def __init__(self, ww3_fname):
        self.fname=ww3_fname
        import numpy as np
        import math
        import datetime
        
        with open(ww3_fname,'r') as file:
            lines = [line.rstrip('\n').split() for line in file]

        nfreq = int(lines[0][3])
        ndir = int(lines[0][4])
        nloc = int(lines[0][5])

        n_spc_per_line = 7
        num=0
        self.freqs=[]
        self.dirs=[]
        for i,line in enumerate(lines[1:]):
            num+=len(line)
            if num<=nfreq:
                for j in line:
                    self.freqs.append(float(j))
                if num == nfreq:
                    n_freq_line=i+1
            elif num<=nfreq+ndir:
                for j in line:
                    self.dirs.append(float(j))
                if num == nfreq+ndir:
                    n_dir_line=i+1-n_freq_line
                    break
            else:
                break

        nhed=n_freq_line+n_dir_line+1
        ntable = int((len(lines)-nhed)/(math.ceil(nfreq*ndir/n_spc_per_line+1)*nloc+1))
        nline_in_table = math.ceil(nfreq*ndir/n_spc_per_line+1)*nloc+1
        nline_in_loc_table = math.ceil(nfreq*ndir/n_spc_per_line)+1

        self.locs = []
        table = lines[nhed:(nhed+nline_in_table)]
        for j in range(nloc):
            line_index1 = 1+j*nline_in_loc_table
            loc_table = table[line_index1:(line_index1+nline_in_loc_table)]
            self.locs.append([float(loc_table[0][3]),float(loc_table[0][2])])
        #######################再用一个大循环提取和写入频谱spc#########################
        self.times=[]
        spcs=[]
        for i in range(ntable): # 每个时间步
            line_index = nhed+i*nline_in_table
            table = lines[line_index:(line_index+nline_in_table)]
            self.times.append(datetime.datetime.strptime(table[0][0]+table[0][1],'%Y%m%d%H%M%S'))
            spc=[]
            for j in range(nloc): # 每个位置
                line_index1 = 1+j*nline_in_loc_table
                loc_table = table[line_index1:(line_index1+nline_in_loc_table)]
                tem=np.array([float(string2) for string1 in loc_table[1:] for string2 in string1]).reshape(ndir,nfreq)
                spc.append(tem.T)
            spcs.append(spc)
        self.spcs=np.array(spcs)#四个维度分别是时间，位置，频率，方向
        
    def spc(self,nstep,nloc):
        return self.spcs[nstep,nloc,:,:]

    def plot_spc(self,nstep,nloc,ax=None,cmap='viridis',arrowdir=45,ymax=None,isave=False,fname=None):
        import matplotlib.pyplot as plt
        import numpy as np
        import math
        import matplotlib
        
        plt.rcParams['font.family'] = 'Times New Roman'
        
        y=self.freqs
        x=self.dirs
        z=self.spcs[nstep,nloc,:,:]
        X,Y=np.meshgrid(x,y)
        # z=z.reshape(z.shape[1],-1).T
        theta_offset=np.pi*90/180
        if ax==None:
            ax=plt.subplot(polar=True,theta_offset=theta_offset)
            ax.set_theta_direction(-1)
        else:
            if not str(type(ax))=="<class 'matplotlib.projections.polar.PolarAxes'>":
                import warnings
                warnings.warn("Given ax is not a polar ax. A new empty ax will be created.",UserWarning)
                ax=plt.subplot(polar=True,theta_offset=theta_offset)
                ax.set_theta_direction(-1)
        # plt.style.use('ggplot')
        if type(cmap)==str:
            cmap=plt.colormaps[cmap]
        if ymax==None:
            ymax=math.ceil(max(Y[z>0.05]*10))/10
        # norm = matplotlib.colors.Normalize(vmin=0.005, vmax=z.max())
        if math.ceil(z.max()*100)/100<0.05:
            levels=np.linspace(0.01,math.ceil(z.max()*100)/100,9)
        else:
            levels=np.linspace(0.05,math.ceil(z.max()*100)/100,9)
        con=ax.contourf(X,Y,z,cmap=cmap,levels=levels)
        ax.set_ylim([0,ymax])
        ax.set_xlabel('directions')
        plt.thetagrids(np.linspace(0,330,12))
        plt.rgrids(np.linspace(0,ymax,3),angle=arrowdir,ha='right',color='black')

        cb=plt.colorbar(con)
        cb.ax.set_ylabel("wave spectrum [m$^2$/Hz/rad]")
        ax.text(np.pi*(arrowdir-5)/180,ymax*0.4,'frequencies [Hz]',ha='center',rotation=-(arrowdir-theta_offset*180/np.pi),color='black')
        if isave:
            if fname==None:
                fname=self.fname[:-3]+f'step{nstep}loc{nloc}.svg'
                plt.savefig(fname,format = 'svg',dpi=200)
                
    def sort_locations(self,seq=[1,2,3,4,5,6,7,9,8,10,11],new_ww3_fname=None):
        import numpy as np
        import math
        
        with open(self.fname,'r') as f:
            lines=f.readlines()
        nfreq = int(lines[0].split()[3])
        ndir = int(lines[0].split()[4])
        nloc = int(lines[0].split()[5])

        n_spc_per_line = 7
        num=0
        for i,line in enumerate(lines[1:100]):
            line=line.split()
            num+=len(line)
            if num<=nfreq:
                if num == nfreq:
                    n_freq_line=i+1
            elif num<=nfreq+ndir:
                if num == nfreq+ndir:
                    n_dir_line=i+1-n_freq_line
                    break
            else:
                break

        nhed=n_freq_line+n_dir_line+1
        ntable = int((len(lines)-nhed)/(math.ceil(nfreq*ndir/n_spc_per_line+1)*nloc+1))
        nline_in_table = math.ceil(nfreq*ndir/n_spc_per_line+1)*nloc+1
        nline_in_loc_table = math.ceil(nfreq*ndir/n_spc_per_line)+1

        new_lines=lines

        for i in range(ntable): # 每个时间步
            line_index = nhed+i*nline_in_table
            table = lines[line_index:(line_index+nline_in_table)]
            table1=table[1:]
            grouped_table = [table1[i:i+101] for i in range(0, len(table1), 101)]
            sorted_groups = [group for _, group in sorted(zip(seq, grouped_table))]
            rearranged_table = [item for sublist in sorted_groups for item in sublist]
            rearranged_table.insert(0,table[0])
            new_lines[line_index:(line_index+nline_in_table)]=rearranged_table
        if new_ww3_fname==None:
            new_ww3_fname=self.fname[:-4]+'_new.spc'
        with open(new_ww3_fname,'w') as f:
            f.writelines(new_lines)
            
    def spc_convert(self):
        # 返回3个ntimes*nlocs的矩阵，分别表示三个变量Hs,Tp,peak_direction（有效波高、峰值周期、峰值方向）。矩阵中的行列表示时间和位置。
        import numpy as np
        frequency_bins = 1/np.array(self.freqs)
        direction_bins = self.dirs
        ntime=len(self.times)
        nloc=len(self.locs)
        
        Hs=np.zeros((ntime,nloc))
        Tp=np.zeros((ntime,nloc))
        peak_direction=np.zeros((ntime,nloc))
        
        for i in range(ntime):
            for j in range(nloc):
                spectrum=self.spc(i,j)
                Hs[i,j] = 4 * np.sqrt(np.trapz(np.trapz(spectrum, direction_bins, axis=1), frequency_bins))
                max_energy_index = np.unravel_index(np.argmax(spectrum, axis=None), spectrum.shape)
                Tp[i,j] = frequency_bins[max_energy_index[0]]
                peak_direction[i,j] = direction_bins[max_energy_index[1]]

                # print("Significant Wave Height (Hs):", Hs)
                # print("Peak Wave Period (Tp):", Tp)
                # print("Peak Wave Direction (Radians):", peak_direction)
        
        return Hs,Tp,peak_direction
    
            
    
def compare_vario_models(samples_xlsx=r"E:\OneDrive\OneDrive - stu.ouc.edu.cn\document\13_我的数据\20231012 block3测量\20231012 block3测量.xlsx",
                        samples_veg_edge_csv=r'X:\GEOSTAT\藨草区域20230105岸线.csv',
                        density_fac = 100,ispherical_gird=True,iplot=True):
    '''
        此脚本用于对比各植被参数的半方差模型拟合情况。
        # 参数含义：
        samples_xlsx 植被样本文件samples_xlsx，格式与r"E:\OneDrive\OneDrive - stu.ouc.edu.cn\document\13_我的数据\20231012 block3测量\20231012 block3测量.xlsx",sheet_name="植被参数"相同
        samples_veg_edge_csv 植被样本区域的岸线文件samples_veg_edge_csv（可以与整个研究区域的岸线相同），格式参考'X:/GEOSTAT/藨草区域20230105岸线.csv'
        density_fac表示植被密度应该乘以多少，才能将单位换算成平方米，例如：采样的样方为0.1m*0.1m，那么fac应当设为100：
            100=1/0.01；如果采样样方为0.2m*0.2m，那么fac=25。
        ispherical_gird表示网格是地理坐标还是笛卡尔坐标。如果半方差模型拟合不好的话，调整line174的scale_length和lag。
        
        eg. :
            generate_plant_rf()
    '''
    import numpy as np
    from matplotlib import pyplot as plt
    from matplotlib.pyplot import get_cmap
    from matplotlib import rcParams
    import matplotlib
    import pandas as pd
    import gstools as gs
    import textwrap
    from scipy.optimize import fmin_cobyla 
    import math
    import my_utils
    
    colors=my_utils.my_colors('gradual_change',n_colors=10,color_set='green2brown')
    rcParams['font.family'] = 'Times New Roman'
    rcParams['font.size'] = 15
    
    # 第一步：读取植被样本
    # 读取文件并生成gs.field类型数据，如果植被样本文件格式不一样，应修改这里的代码
    df = pd.read_excel(samples_xlsx,sheet_name="Sheet1")
    data_array = np.array(df.set_index('编号'))
    lat=data_array[:,1]
    lon=data_array[:,2]
    x=data_array[:,3]
    y=data_array[:,4]
    ds=data_array[:,5]*density_fac/100
    ht=data_array[:,6]
    ht[0]=np.mean(ht[1::])
    dm=data_array[:,7]
    scale_lon=x.max()-x.min()
    scale_lat=y.max()-y.min()
    # if ispherical_gird:
        # coord_unit='[°]'
        # scale_length=0.0011
        # lag=0.0001
        
    # else:
        # lon=x
        # lat=y
        # coord_unit='[m]'
        # scale_length=110
        # lag=12.35
        
    if not ispherical_gird:
        lon=x
        lat=y
        coord_unit='[m]'

    else:
        coord_unit='[°]'
        
    scale_length=[120,120,120]
    lag=[9.9,12.2,13.5]
        
    # 把距离标准化到0~scale范围内：
    def std_lonlat(lonlat,reverse=False,ori_lonlat=None,scale=500):
        if not reverse:
            mmin=lonlat.min()
            mmax=lonlat.max()
            lonlat1=scale*(lonlat-mmin)/(mmax-mmin)
        else:
            mmin=ori_lonlat.min()
            mmax=ori_lonlat.max()
            lonlat1=lonlat/scale*(mmax-mmin)+mmin
        return lonlat1
    ori_lon=lon
    ori_lat=lat
    
    # 第二步：拟合趋势（这里拟合的是一维趋势）
    # 读取岸线
    
    land_bound=pd.read_csv(samples_veg_edge_csv,header=None).values
    
    #计算样本点和海岸线的最短距离
    target_points=np.array([lon,lat]).T
    sampledis_to_csl=np.array([my_utils.calculate_min_dis(target_points[i,:],land_bound) for i in range(target_points.shape[0])])
    def plt_linear_trend(x,y,ylabel,xlabel='Distance to vegetation edge '+coord_unit,colors=[],ax=[],iplot=True):
        from scipy.stats import t
        if len(colors)==0:
            colors=my_utils.my_colors('contrast',n_colors=10)
        # 设置字体和字号
        # rcParams['font.family'] = 'Times New Roman'
        # plt.rcParams['font.family'] = 'Times New Roman'
        # plt.rcParams['font.size'] = 15
        # 使用 polyfit() 函数进行线性拟合
        coefficients = np.polyfit(x, y, 1)
        slope = coefficients[0]
        intercept = coefficients[1]
        # 计算拟合值
        predicted = slope * x + intercept
        # 计算R2
        r_squared = 1 - (np.sum((y - predicted)**2) / np.sum((y - np.mean(y))**2))
        
        # 得到斜率和截距的标准误差
        n = len(x)
        std_error_slope = np.sqrt(np.sum((y - predicted)**2) / (n - 2)) / np.sqrt(np.sum((x - np.mean(x))**2))
        std_error_intercept = std_error_slope * np.sqrt(np.sum(x**2) / n)
        # 计算置信区间（以95%置信水平为例）
        alpha = 0.05  # 置信水平为95%
        t_critical = t.ppf(1 - alpha/2, n - 2)
        slope_lower, slope_upper = slope - t_critical * std_error_slope, slope + t_critical * std_error_slope
        intercept_lower, intercept_upper = intercept - t_critical * std_error_intercept, intercept + t_critical * std_error_intercept
        if iplot:
            print(f"拟合的R平方值：{r_squared}")
            print(f"斜率置信区间：({slope_lower}, {slope_upper})")
            print(f"截距置信区间：({intercept_lower}, {intercept_upper})")
        # # 绘制原始数据和拟合直线以及置信区间
        # if iplot:
            if ax==[]:
                ax = plt.gca()
            ax.scatter(x, y, label='samples',marker='s',color=colors[0])
            ax.plot(x, slope * x + intercept, color=colors[3], label='fitted')
            xx = np.linspace(0.95*min(x), 1.05*max(x), 10)
            ax.fill_between(xx, (slope_lower * xx + intercept_lower), (slope_upper * xx + intercept_upper), color=colors[0], alpha=0.2, label='95% CI')
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.legend()
            ax.text(0.13,0.82,f'R$^2$ :  {r_squared:5f}',transform=ax.transAxes)
            # plt.savefig(ylabel+'.png',dpi=200)
            # plt.show()
        return r_squared, slope, intercept
    slope = np.zeros(3)
    intercept = np.zeros(3)
    r_squared = np.zeros(3)
    if iplot:
        fig, axes = plt.subplots(1, 3, figsize=[18,5.5],tight_layout=True)
    else:
        axes=[[],[],[]]

    r_squared[0], slope[0],intercept[0] = plt_linear_trend(sampledis_to_csl,ds,"density [stems/0.01m$^2$]",ax=axes[0],iplot=iplot)#,colors=colors)
    r_squared[1], slope[1],intercept[1] = plt_linear_trend(sampledis_to_csl,ht,"height [cm]",ax=axes[1],iplot=iplot)#,colors=colors)
    r_squared[2], slope[2],intercept[2] = plt_linear_trend(sampledis_to_csl,dm,"diameter [mm]",ax=axes[2],iplot=iplot)#,colors=colors)
    

    if iplot:
        plt.savefig('3par.png',dpi=200)
        plt.show()

    # 第三步：拟合半方差模型，以°为距离单位
    
    slon=std_lonlat(lon,scale=scale_lon)
    slat=std_lonlat(lat,scale=scale_lat)
    
    def trend(x,y,coast_line,slope,intercept):
        target_points=np.array([x,y]).T
        dis=np.array([my_utils.calculate_min_dis(target_points[i,:],coast_line) for i in range(target_points.shape[0])])
        trend1=dis*slope+intercept
        return trend1
    for i,var in enumerate(['ds','ht','dm']):
        if r_squared[i]>0.1:
            exec(f"{var}_resi={var}-trend(lon,lat,land_bound,slope[{i}],intercept[{i}])")
        else:
            exec(f"{var}_resi={var}-np.mean({var})")
    # ds_resi=ds-np.mean(ds)
    # ht_resi=ht-trend(lon,lat,land_bound,slope[1],intercept[1])
    # dm_resi=dm-np.mean(dm)

    ################# 循环定义变量##################
    for i,var_name in enumerate(['ds','ht','dm','ds_resi','ht_resi','dm_resi']):
        exec("%s_field=gs.field.Field(dim=2)"%var_name)
        exec("%s_field(pos=(x,y),field=%s)"%(var_name,var_name))

    def test_vario_model(xy,field,max_lag,bin_wids,models,ishow=True,ax=[],colors=my_utils.my_colors('contrast',n_colors=10),linestyles=['-','--'],return_ax=False,return_bin=False,return_emp_var=False):
        
        scores = {}
        paras = {}
        pcovs = {}
        fit_models = {}
        
        bins = np.arange(0, max_lag, bin_wids)
        bin_center, vario = gs.vario_estimate(
            *(xy, field, bins),
            mesh_type="unstructured"
        )
        if ishow:
            if ax==[]:
                ax = plt.gca()
            # sca=ax.scatter(bin_center, vario, color=colors[9], label="empirical",marker='s')

            
            
        for i,model in enumerate(models):
            fit_models[model] = models[model](dim=2)
            para, pcov, r2 = fit_models[model].fit_variogram(bin_center, vario, return_r2=True)
            scores[model] = r2
            paras[model] = para
            pcovs[model] = pcov
            if ishow:
                fit_models[model].plot(func='vario_yadrenko',x_max=max_lag, ax=ax,color=colors[i],linestyle=linestyles[i//(len(models)//2)])
                # ax.legend(loc="lower right")
        ranking = sorted(scores.items(), key=lambda item: item[1], reverse=True)
        if ishow:
            print("RANKING by Pseudo-r2 score")
            for i, (model, score) in enumerate(ranking, 1):
                print(f"{i:>6}. {model:>15}: {score:.5}")
        best_model_type = ranking[0][0]
        best_model = fit_models[best_model_type]

        return paras[best_model_type], pcovs[best_model_type], scores, best_model,ax,bin_center,vario

    ######绘制拟合半方差模型的图像
    models = {
        "Gaussian": gs.Gaussian,
        "Exponential": gs.Exponential,
        "Matern": gs.Matern,
        "Stable": gs.Stable,
        "Rational": gs.Rational,
        "Circular": gs.Circular,
        "Spherical": gs.Spherical,
        "SuperSpherical": gs.SuperSpherical,
        #"JBessel": gs.JBessel,
    }
    # variables=['ds','ht','dm','ds_resi','ht_resi','dm_resi']
    variables=['ds_resi','ht_resi','dm_resi']
    units=[' [stems/m$^2$]',' [cm]',' [mm]',' [stems/0.01m$^2$]',' [cm]',' [mm]']
    var_long_name=['density','height','diamter','density','height','diamter']

    # paras={}
    # pcovs={}
    r2s={}
    bst_mds={}
    if iplot:
        rcParams['font.size'] = 7
        fig, axes = plt.subplots(1, 3, figsize=[12, 3],tight_layout=True)

    for i,var in enumerate(variables):
        _,_,r2s[var],bst_mds[var] ,_,_,_= test_vario_model((slon,slat),eval(var+'_field'),scale_length[i],lag[i],models,ax=axes[np.mod(i,3)],ishow=iplot)#i//3,
        if iplot:
            axes[np.mod(i,3)].set_title('fitting plant '+var_long_name[i]+' to an isotropic model')
            axes[np.mod(i,3)].set_xlabel('distance '+'[m]') #coord_unit
            axes[np.mod(i,3)].set_ylabel('semivariance')
            txt='        R$^2$: \n'
            for model,r2 in r2s[var].items():
                txt=txt+f'{model:>15}: {r2:<.5}\n'
            axes[np.mod(i,3)].text(0.3,0,txt,transform=axes[np.mod(i,3)].transAxes, ha='left')
            if i==len(variables)-1:
                axes[np.mod(i,3)].legend(loc="lower right")
    if iplot:
        plt.savefig('compare_vario_models.png',dpi=200)