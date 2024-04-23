import xarray as xr
import numpy as np
import pandas as pd
import datetime

@np.vectorize
def in_range(value, start, stop):
    return start <= value <= stop
    
def my_colors(style, n_colors=5,color_set='green2brown'):
    '''Use my_utils.try_colors(<color list of string (hex index)>)to see what colors you have got.'''
    # warning: the default all_colors just have 10 colors, n_colors should be less than or equal to 10,
    import matplotlib.pyplot as plt
    import numpy as np
    if color_set=='green2brown':
        all_colors=['#017987','#2b9ca1','#7bb2b0','#a4c6bb','#f2dabe','#c7a37e','#7d674b','#332d15'] # '#dbe9e3','#fff7ee'太浅了，如果需要的话再插在'#a4c6bb','#f2dabe'之间
    elif color_set=='blue2red':
        all_colors=['#4257b0','#9be4d9','#fff1a6','#d64951','#9b9cbb']
    elif color_set=='classic':
        all_colors=['#202781','#3451a5','#396db5','#46c7ec','#82ca99','#bed735','#fbce09','#ef5e1b','#ec171c','#7e0b0d']
    elif color_set=='contrast':
        all_colors=['#f2c909','#fd625e','#344447','#00b9ac']
    else:
        raise ValueError(f"There is no colorset name {color_set}, please choose from 'green2brown', 'blue2red', 'classic', 'contrast'.")

    if style == 'contrast':
        all_colors=np.array(all_colors)
        ind = np.linspace(0,len(all_colors)-1,n_colors,dtype='int32')
        ind1 = np.append(np.append(ind[::3],ind[1::3]),ind[2::3])
        out_colors = all_colors[ind1]
    elif style == 'gradual_change':
        all_colors=np.array(all_colors)
        # ind = np.linspace(0,len(all_colors),n_colors,dtype='int32')
        if n_colors > len(all_colors)/2:
            out_colors=all_colors[0:n_colors]
        else:
            out_colors = all_colors[0::2][0:n_colors]
    elif style == 'colormap':
        import matplotlib.colors as colors
        out_colors = colors.LinearSegmentedColormap.from_list("brw", all_colors, N=256)
    else:
        raise ValueError("Invalid style. Choose from 'contrast', 'gradual_change', or 'colormap'.")
    return out_colors

def try_colors(colors=['#017987','#2b9ca1','#7bb2b0','#a4c6bb','#dbe9e3','#fff7ee','#f2dabe','#c7a37e','#7d674b','#332d15']):
    import matplotlib.pyplot as plt
    x = [1, 2, 3, 4, 5]
    y = [10, 15, 7, 12, 8]
    for i,color in enumerate(colors):
        plt.plot(x, [j-2*i for j in y], color=color, label=color)
    plt.legend(bbox_to_anchor=(1.1,0.5),loc=6, borderaxespad=0)
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.title('Try your own colors')
    plt.show()
    # usage example:
    # my_utils.try_colors(my_utils.my_colors('gradual_change',n_colors=10))
    # my_utils.try_colors(my_utils.my_colors('contrast',n_colors=10))


def calculate_min_dis(target_point,lines):
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
        return distance(px, py, perpendicular[0], perpendicular[1])
    # 初始化最小距离为一个很大的数
    min_distance = float('inf')
    # 将点集合转换为线段集合
    segments=points_to_segments(lines)
    # 遍历每个线段，计算最短距离
    for segment in segments:
        x1, y1, x2, y2 = segment
        distance = point_to_segment_distance(target_point[0], target_point[1], x1, y1, x2, y2)
        min_distance = min(min_distance, distance)
    return min_distance
    
def polygon_sort(pp): 
    import math
    '''
    按照多边形顶点顺序对数组排序
    '''
    cent=(sum([p[0] for p in pp])/len(pp),sum([p[1] for p in pp])/len(pp))
    pp.sort(key=lambda p: math.atan2(p[1]-cent[1],p[0]-cent[0]))
    return pp

def get_new_block(block1,block2):
    '''
    删除两个n*2数组中重复的元素（两个列表中都要删除），合并两个数组，并且按照多边形顶点顺序对数组排序
    '''
    block1_list = block1.tolist()
    block2_list = block2.tolist()
    repeat_element = []
    [repeat_element.append(x) for x in block1_list if x in block2_list]
    res = []
    [res.append(x) for x in block1_list if x not in repeat_element]
    [res.append(x) for x in block2_list if x not in repeat_element]
    # res = np.array(polygon_sort(res))
    return np.array(res)
    
def create_proj_dir(proj_name, nnest=1, homedir=r'E:\d3d_cases', subfolders=['/1_swn','/2_wind','/3_grid','/4_bound','/5_coupled_flow','/6_attribute_files','/7_tmp']):
    '''
        this function is used to create project folder tree.
        the main folder name is the same as proj_name; 
        nnest is the number of nested layers. nnest=1 meas that there is no nesting (just one layer).
        if nnest>1,there will be subfolder of each layer;
        and then in each layer, subfolders would be created by the parameter 'subfolders'.
        
        example:create_proj_dir('test',nnest=2)
    '''
    # writing & debugging
    # proj_name='ttt'
    # nnest=1
    # homedir=r'E:\d3d_cases'
    # subfolders=['final_control','wind','grid','bound','attribute_files','obs']
    
    import os
    # initialize
    proj_dir=rf'{homedir}\{proj_name}'
    if os.path.exists(proj_dir):
        print("文件夹已存在，是否覆盖？")
        while True:
            user_input = input("请输入'y'或者'n'：")
            if user_input == 'y':
                import shutil
                shutil.rmtree(proj_dir)
                break
            elif user_input == 'n':
                print('程序退出')
                return
            else:
                print('无效的输入！请重新输入。')
    
    # creat folders
    os.mkdir(proj_dir)
    if nnest==1:
        for folder in subfolders:
            os.mkdir(rf'{proj_dir}\{folder}')
    else:
        os.mkdir(rf'{proj_dir}\outer')
        for folder in subfolders:
            os.mkdir(rf'{proj_dir}\outer\{folder}')
        if nnest==2:
            os.mkdir(rf'{proj_dir}\inner')
            for folder in subfolders:
                os.mkdir(rf'{proj_dir}\inner\{folder}')
        else:
            for i in range(nnest-1):
                j=i+1
                os.mkdir(rf'{proj_dir}\inner{j}')
                for folder in subfolders:
                    os.mkdir(rf'{proj_dir}\inner{j}\{folder}')
                    
def format_str(data,n_data_per_line,formater=' 1.7E',nspace=1):
    stri=''
    ndata=len(data)
    a="f'{x:"+formater
    b="}'"
    c=a+b
    if ndata%n_data_per_line==0:
        for i in range(ndata//n_data_per_line):
            for j in range(n_data_per_line):
                x=data[i*n_data_per_line+j]
                stri=stri+nspace*" "+(eval(c))
            stri=stri+'\n'
    else:
        for i in range(ndata//n_data_per_line):
            for j in range(n_data_per_line):
                x=data[i*n_data_per_line+j]
                stri=stri+nspace*" "+(eval(c))
            stri=stri+'\n'
        for j in range(ndata%n_data_per_line):
            x=data[(ndata//n_data_per_line)*n_data_per_line+j]
            stri=stri+nspace*" "+(eval(c))
        stri=stri+'\n'
    return stri
    
def write_formated_file(array,fname,formater=' 1.7E',max_num_in_line=12):
    a="f'{elem:"+formater
    b="}  '"
    c=a+b

    with open(fname, "w") as f:
        for row in array:
            count = 0
            for elem in row:

                f.write(eval(c))
                count += 1
                if count == max_num_in_line: # 在这里改每行最多几个数字
                    f.write('\n')
                    count = 0
            if count > 0:
                f.write('\n')

def replace_strs_in_file(filenm,ori_strs,new_strs,new_fnm=[]):
    if type(ori_strs)==str:
        ori_strs=[ori_strs]
        new_strs=[new_strs]
    if new_fnm==[]:
        tmp=filenm.split('.')[:-1]
        tmp.append('_new.')
        tmp.append(filenm.split('.')[-1])
        new_fnm=''.join(tmp)
    with open(filenm, 'r') as f:
        lines=f.read()
    for i,ori_str in enumerate(ori_strs):
        lines=lines.replace(ori_str,new_strs[i]) 
    with open(new_fnm, 'w') as f:
        f.write(lines)