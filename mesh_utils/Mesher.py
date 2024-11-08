# -*- coding: UTF-8 -*-
'''
根据散点的数据生成三维面，包括多段线和散点

'''
import numpy as np
import matplotlib.tri as mtri
from read import F
from tool import binary_to_decimal

# 需要从UI输入的参数
# 设定计算场地范围，这个范围必须在的散点范围内，不能在散点范围外
minx = 0.0
maxx = 1000.0
miny = 0.0
maxy = 1500.0
numx = 9 #x方向上需要剖分几个网格
numy = 8 #y方向上需要剖分几个网格
all_sizez = [3,2,3,3,3] #每个材料层在z方向上剖分几个网格
path = '/home/momo/geos/mesh_utils/' # 文件路径
file = 'points-large.in' #文件名


# 计算
NN = numx * numy
count = 0
# 读取散点文件
In = F(path, file)
In.Read()
matnum = In.matnum
# 按照输入的参数，在xy平面上形成meshgrid
Nx =  np.linspace(minx, maxx, numx)
Ny = np.linspace(miny, maxy, numy)
gridx, gridy = np.meshgrid(Nx,Ny)
# 计算最表面的一层
surface = In.mat[0]
x0 = surface['X']
y0 = surface['Y']
z0 = surface['Z']
triang0 = mtri.Triangulation(x0, y0)
interp_z0 = mtri.LinearTriInterpolator(triang0, z0)
zi0 = interp_z0(gridx,gridy)
Nnode = 0
NnodeinTP = 0
Ntp = 0
NdType = []
with open('node.out','w') as f:#初始化nodeout文件，存储网格结点信息，并且写入最表层的node
    f.write('#NodeNum----x----y----z'+"\n")
    for j in range(numy):
        for i in range(numx):
            Nnode = Nnode + 1
            xx = gridx[j][i]
            yy = gridy[j][i]
            f.write(str(Nnode) + "\t")
            f.write("{:.{}e}".format(xx,3) + "\t")
            f.write("{:.{}e}".format(yy,3) + "\t")
            f.write("{:.{}e}".format(zi0[j][i],3) + "\t")
            f.write("0"+"\t")
            typebinary = "000011"
            if xx == minx:
                typebinary = '10' + typebinary[2:]
            if xx == maxx:
                typebinary = '11' + typebinary[2:]
            if yy == miny:
                typebinary = typebinary[:2] + '10' + typebinary[4:]
            if yy == maxy:
                typebinary = typebinary[:2] + '11' + typebinary[4:]
            typed = str(binary_to_decimal(typebinary))
            NdType.append(typebinary)
            f.write(typed+"\n")
            
f.close()
with open('mesh.out','w') as m:#初始化meshout文件，存储网格拓扑信息，包括材料编号
    m.write("#No"+"\t"+"MAT"+"\t"
            +"node1"+"\t"
            +"node2"+"\t"
            +"node3"+"\t"
            +"node4"+"\t"
            +"node5"+"\t"
            +"node6"+"\t"
            +"node7"+"\t"
            +"node8"+"\t"+"\n")
m.close()
with open('notion.out','w') as n:
    print('notion created')
n.close()
all_data = In.mat
x0 = all_data[0]['X']
y0 = all_data[0]['Y']
z0 = all_data[0]['Z']
triang0 = mtri.Triangulation(x0, y0)
interp_z0 = mtri.LinearTriInterpolator(triang0, z0)
zi0 = interp_z0(gridx,gridy)
all_data.pop(0)
for index,lower in enumerate(all_data):
    x = lower['X']
    y = lower['Y']
    z = lower['Z']
    Nmat = lower['mat']
    sizez = all_sizez[Nmat - 1]
    triang = mtri.Triangulation(x, y)
    interp_z = mtri.LinearTriInterpolator(triang, z)
    zi = interp_z(gridx, gridy)
    with open('node.out','a') as f ,open('mesh.out','a') as m,open('notion.out','a') as n:
        for k in range(sizez):
            for j in range(numy):
                for i in range(numx):
                    #结点
                    xx = gridx[j][i]
                    yy = gridy[j][i]
                    Nnode = Nnode + 1
                    f.write(str(Nnode) + "\t")
                    f.write("{:.{}e}".format(xx,3) + "\t")
                    f.write("{:.{}e}".format(yy,3) + "\t")
                    zi1 = zi[j][i]
                    zi2 = zi0[j][i]
                    zz = zi2-(k+1)*(zi2-zi1)/sizez
                    f.write("{:.{}e}".format(zz,3) + "\t")
                    f.write(str(Nmat)+"\t")
                    typebinary = "000000"
                    if xx == minx:
                        typebinary = '10' + typebinary[2:]
                    if xx == maxx:
                        typebinary = '11' + typebinary[2:]
                    if yy == miny:
                        typebinary = typebinary[:2] + '10' + typebinary[4:]
                    if yy == maxy:
                        typebinary = typebinary[:2] + '11' + typebinary[4:]
                    if k == sizez-1 and index == len(all_data)-1:
                        typebinary = typebinary[:4] + '10'
                    NdType.append(typebinary)
                    typed = str(binary_to_decimal(typebinary))
                    f.write(typed)
                    f.write('\n')
                    n.write(typed)
                    n.write('\n')
                    #拓扑
                    NnodeinTP = NnodeinTP +1
                    if i == numx - 1 or j == numy - 1:
                        continue
                    Ntp = Ntp + 1
                    m.write(str(Ntp)+"\t") # 单元编号
                    m.write(str(Nmat)+"\t") # 材料编号
                    m.write(str(NnodeinTP)+"\t") # 结点1
                    m.write(str(NnodeinTP + 1)+"\t") # 结点2
                    m.write(str(NnodeinTP + numx +1)+"\t") # 结点3
                    m.write(str(NnodeinTP + numx)+"\t") # 结点4
                    m.write(str(NnodeinTP + NN)+"\t") # 结点5
                    m.write(str(NnodeinTP + NN + 1)+"\t") # 结点6
                    m.write(str(NnodeinTP + NN + numx + 1)+"\t") # 结点7
                    m.write(str(NnodeinTP + NN + numx)+"\n") # 结点8
                    
    x0 = x
    y0 = y
    z0 = z
    triang0 = mtri.Triangulation(x0, y0)
    interp_z0 = mtri.LinearTriInterpolator(triang0, z0)
    zi0 = interp_z0(gridx,gridy)

tecplot = True
if tecplot:
    print('making tecplt file')
    content1 = []
    content2 = []
    with open('node.out','r') as f:
        for line in f:
            if line[0] == '#': # 跳过标识行
                continue
            line = line[0:-1]
            line = line.split()
            line.pop(0)
            content1.append(' '.join(line))
    with open('mesh.out','r') as f:
        for line in f:
            if line[0] == '#': # 跳过标识行
                continue
            line = line[0:-1]
            line = line.split()
            line.pop(0)
            line.pop(0)
            content2.append( ' '.join(line))
    ttnd = len(content1)
    tttp = len(content2)
    with open('tecplt.dat','w') as t:
        t.write('TITLE     = "Tecplot"')
        t.write('\n')
        t.write('VARIABLES = "X(m)"')
        t.write('\n')
        t.write('"Y(m)"')
        t.write('\n')
        t.write('"Z(m)"')
        t.write('\n')
        t.write('"MAT"')
        t.write('\n')
        t.write('"BC"')
        t.write('\n')
        t.write('ZONE T="MI= 1  MAT=  1"')
        t.write('\n')
        t.write('STRANDID=0, SOLUTIONTIME=0')
        t.write('\n')
        t.write(f'Nodes= {ttnd}, Elements=  {tttp}, ZONETYPE=FEBrick')
        t.write('\n')
        t.write('DATAPACKING=POINT')
        t.write('\n')
        t.write('DT=(SINGLE SINGLE SINGLE )')
        t.write('\n')
        for item in content1:
            t.write(item)
            t.write('\n')
        for item in content2:
            t.write(item)
            t.write('\n')

print('tecplot file is ready')

vtk = True
if vtk:
    temp_dic = {}
    typecode = 80
    print('making vtk file')
    content1 = []
    content2 = []
    content3 = []
    with open('node.out','r') as f:
        for line in f:
            if line[0] == '#': # 跳过标识行
                continue
            line = line[0:-1]
            line = line.split()
            line.pop(0)
            line.pop(-1) 
            line.pop(-1) 
            content1.append(' '.join(line))
    with open('mesh.out','r') as f:
        for line in f:
            if line[0] == '#': # 跳过标识行
                continue
            line = line[0:-1]
            line = line.split()
            line.pop(0)
            line.pop(0)
            incremented_line = [str(int(item) - 1) for item in line]
            content2.append( ' '.join(incremented_line))
    ttnd = len(content1)
    tttp = len(content2)
    with open('mesh.vtk','w') as t:
        t.write('# vtk DataFile Version 2.0\n')
        t.write('Created by MOMOMesher\n')
        t.write('ASCII\n')
        t.write('DATASET UNSTRUCTURED_GRID\n')
        t.write(f'POINTS {ttnd} double\n')
        for item in content1:
            t.write(item)
            t.write('\n')
        t.write(f'CELLS {tttp} \t {tttp*9} \n')
        for item in content2:
            t.write('8\t')
            t.write(item)
            t.write('\n')
        t.write(f'CELL_TYPES {tttp} \n')
        for item in content2:
            t.write('12\n')
        t.write(f'CELL_DATA {tttp} \n')
        t.write('SCALARS CellEntityIds int 1 \n')
        t.write('LOOKUP_TABLE defult \n')
        for item in content2:
            item = item.split()
            typeII = ''
            for nd in item:
                nd = int(nd)
                typeII=typeII + NdType[nd]
            if typeII in temp_dic.keys():
                t.write(f'{temp_dic[typeII]}\n')
            else:
                temp_dic[typeII] = typecode
                t.write(f'{typecode}\n')
                typecode = typecode + 1


print('end')
