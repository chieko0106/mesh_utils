# -*- coding: UTF-8 -*-
class F:
    mat = []
    matnum = 0
    data = {}
    count = 0
    def __init__(self,path,file):
        self.__path = path
        self.__file = file

    def Read(self):
        __f = open(self.__path + self.__file)
        for line in __f:
            if line == '%END':
                F.mat.append(F.data)
                F.matnum = F.matnum - 1
                print('%END detected, end read infile')
                break
            if line[0] == '#': # 跳过标识行
                continue
            if line[0] == '%': # 新材料标识
                if F.matnum != 0:
                    F.mat.append(F.data)
                F.data = {'X':[],'Y':[],'Z':[],'count':0,'mat':F.matnum}
                F.matnum = F.matnum + 1
                continue
            if len(line) < 1:# 跳过空行
                print('warning: in infile, empty line detected')
                continue
            line = line[0:-1]
            line = line.split()
            F.data['X'].append(float(line[0]))
            F.data['Y'].append(float(line[1]))
            F.data['Z'].append(float(line[2]))
            F.data['count'] += 1
    def maskRead(self,X,Y,Z):
        xx = F.data['X']
        yy = F.data['Y']
        zz = F.data['Z']
        x = []
        y = []
        z = []
        F.count = 0
        for i in range(len(xx)):
            if xx[i] > X[0] and xx[i] < X[1] and yy[i] > Y[0] and yy[i] < Y[1]and zz [i] > Z[0] and zz[i] < Z[1]:
                x.append(xx[i])
                y.append(yy[i])
                z.append(zz[i])
                F.count += 1
        F.data['X'] = x
        F.data['Y'] = y
        F.data['Z'] = z
