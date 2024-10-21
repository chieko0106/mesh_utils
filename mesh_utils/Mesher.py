# -*- coding: UTF-8 -*-
'''
根据散点的数据生成三维面，包括多段线和散点

'''
import numpy as np
import matplotlib.tri as mtri
from read import F

# 需要从UI输入的参数
# 设定计算场地范围，这个范围必须在的散点范围内，不能在散点范围外
minx = 0.0
maxx = 1000.0
miny = 0.0
maxy = 1500.0
numx = 100 #x方向上需要剖分几个网格
numy = 150 #y方向上需要剖分几个网格
all_sizez = [10,5,10,10,10] #每个材料层在z方向上剖分几个网格
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
with open('node.out','w') as f:#初始化nodeout文件，存储网格结点信息，并且写入最表层的node
    f.write('#NodeNum----x----y----z'+"\n")
    for j in range(numy):
        for i in range(numx):
            Nnode = Nnode + 1
            f.write(str(Nnode) + "\t")
            f.write("{:.{}e}".format(gridx[j][i],3) + "\t")
            f.write("{:.{}e}".format(gridy[j][i],3) + "\t")
            f.write("{:.{}e}".format(zi0[j][i],3) + "\t")
            f.write("0"+"\n")
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
surface = True
for lower in all_data:
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
                    f.write(str(Nmat)+"\n")

                    if xx == minx:
                        n.write('00')
                        n.write('\n')
                    elif xx == maxx:
                        n.write('01')
                        n.write('\n')
                    elif yy == miny:
                        n.write('10')
                        n.write('\n')
                    elif yy == maxy:
                        n.write('11')
                        n.write('\n') 
                    elif surface and k == 0:
                        n.write('33') # 最表面点标记为33
                        n.write('\n') 
                    else:
                        n.write('99') #内部点标记为99
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
                    
    surface = False
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
        t.write('"MAT(m)"')
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
print('end')
