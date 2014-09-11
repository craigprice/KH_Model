#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import math
from decimal import *
from os import listdir

PI = 3.14159265358979323846;


#Start Main Function
#Get list of all sim files
#pathToSimFiles ='/Volumes/Toshi/dataFromSimulations/honV9/hep2/dataFiles'
pathToSimFiles ='/Volumes/Toshi/dataFromSimulations/copyPlot'
simFiles = [ f for f in listdir(pathToSimFiles) if f.find('.txt') > 0 ]

#For each Sim File extract and make a plot of the FT
for f in simFiles:
    print f
    opName = f[f.find('OP_')+3:-4]
    pieces = f.split('_')
    graphDescription = ''
    
    #Extract sim parameters
    numSpaces = 7
    for p in range(1, len(pieces)-3, 2):
        tempStr = pieces[p+1]
        tempStr = tempStr.rjust(numSpaces - len(pieces[p])) + '\n'
        graphDescription = graphDescription + pieces[p] + ' = ' + tempStr
    #print graphDescription

    print "Open File to get basic characteristics"
    myfile = open(str(pathToSimFiles + '/' + f),'r')
    numConfigsDone = 0
    cellsA = 0
    cellsB = 0
    cellsC = 0
    KbT = 0
    for line in myfile:
        if line.find('Cells-A: ') == 0:
            cellsA = int(line[len('Cells-A: '):])
        if line.find('KbT: ') == 0:
            KbT = float(line[len('KbT: '):])
        if line.find('Cells-B: ') == 0:
            cellsB = int(line[len('Cells-B: '):])
        if line.find('Cells-C: ') == 0:
            cellsC = int(line[len('Cells-C: '):])
        if line.find('numConfigsDone: ') == 0:
            numConfigsDone = int(line[len('numConfigsDone: '):])
    myfile.close()
    if cellsA < 80:
        pass
        #continue
    if KbT > 0.01:
        continue
    #cellsA = cellsA/3#edit
    #cellsB = cellsB/3#edit
    avector_1x = 1*1;
    avector_2x = 1*0.5;
    avector_1y = 1*0;
    avector_2y = 1*math.sqrt(3.0) / 2.0; 
    delta_x = 0
    delta_y = 1.0/math.sqrt(3.0)
    
    xSpinDirection_x = ((-1)*avector_2x) + delta_x
    xSpinDirection_y = ((-1)*avector_2y) + delta_y
    ySpinDirection_x = ((-1)*avector_2x) + avector_1x + delta_x
    ySpinDirection_y = ((-1)*avector_2y) + avector_1y + delta_y
    zSpinDirection_x = delta_x
    zSpinDirection_y = delta_y
    
    rvector_x = []
    rvector_y = []
    spin_x_u = []
    spin_x_v = []
    spin_y_u = []
    spin_y_v = []
    spin_z_u = []
    spin_z_v = []
    x2 = []
    y2 = []
    u2 = []
    v2 = []
    #print str(numConfigsDone) + " " + str(cellsA) + " " + str(cellsB) + " " + str(cellsC) + " " + str(numFTSumsDone)


    print "Extract components data"
    print "sites in X direction, sites in Y direction, sites in Z direction, sublattice, X component of Spin, Y component of Spin, Z component of Spin"
    myfile = open(str(pathToSimFiles + '/' + f),'r')
    for line in myfile:
        if line.find('Components of Spin ') < 0:
            continue
        component = line[len('Components of Spin '):]
        #print line[:-1]
        xPos = int(component[:component.find(","):])
        component = component[1 + component.find(","):]
        yPos = int(component[:component.find(",")])
        component = component[1 + component.find(","):]
        zPos = int(component[:component.find(",")])
        component = component[1 + component.find(","):]
        sPos = int(component[:component.find(" ")])
        component = component[1 + component.find(" "):]
        spin_x = float(component[:component.find(" ")])
        component = component[1 + component.find(" "):]
        spin_y = float(component[:component.find(" ")])
        component = component[1 + component.find(" "):]
        spin_z = float(component[:component.find(" ")])
        #if ((xPos * avector_1x + yPos * avector_2x + sPos * delta_x) / cellsA) < 1.2:#edit
        #    continue;
        #if ((xPos * avector_1y + yPos * avector_2y + sPos * delta_y) / cellsB) < 0.6:#edit
        #    continue;
            
        rvector_x.append((xPos * avector_1x + yPos * avector_2x + sPos * delta_x) / cellsA)
        rvector_y.append((xPos * avector_1y + yPos * avector_2y + sPos * delta_y) / cellsB)
        spin_x_u.append(spin_x * xSpinDirection_x)
        spin_x_v.append(spin_x * xSpinDirection_y)
        spin_y_u.append(spin_y * ySpinDirection_x)
        spin_y_v.append(spin_y * ySpinDirection_y)
        spin_z_u.append(spin_z * zSpinDirection_x)
        spin_z_v.append(spin_z * zSpinDirection_y)
        print xPos, yPos, zPos, sPos, spin_x, spin_y, spin_z
    
    print "Plot"
    
    #w = 16
    #h = 10
    w = 10
    h = 10
    plt.figure(figsize=(w,h))
    plt.axis((0,w/10.0,0,h/10.0))
    #plt.axis((1.2,1.6,0.6,1.0))#edit
    Q = plt.quiver(rvector_x,rvector_y,spin_x_u,spin_x_v,color='red',pivot='middle',headaxislength=1.1,headlength=1.1,scale=60,width=w/10000.0,headwidth=2)
    Q = plt.quiver(rvector_x,rvector_y,spin_y_u,spin_y_v,color='blue',pivot='middle',headaxislength=1.1,headlength=1.1,scale=60,width=w/10000.0,headwidth=2)
    Q = plt.quiver(rvector_x,rvector_y,spin_z_u,spin_z_v,color='green',pivot='middle',headaxislength=1.1,headlength=1.1,scale=60,width=w/10000.0,headwidth=2)
    #Q = plt.quiver(rvector_x,rvector_y,spin_x_u,spin_x_v,color='red',pivot='middle',headaxislength=1.1,headlength=1.1,scale=8,width=w/1000.0,headwidth=2)#edit
    #Q = plt.quiver(rvector_x,rvector_y,spin_y_u,spin_y_v,color='blue',pivot='middle',headaxislength=1.1,headlength=1.1,scale=8,width=w/1000.0,headwidth=2)#edit
    #Q = plt.quiver(rvector_x,rvector_y,spin_z_u,spin_z_v,color='green',pivot='middle',headaxislength=1.1,headlength=1.1,scale=8,width=w/1000.0,headwidth=2)#edit
    plt.tick_params(\
    axis='both',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    right='off',         # ticks along the top edge are off
    left='off',         # ticks along the top edge are off
    labelbottom='off',
    labelleft='off')#edit    
    plt.savefig('figuresSnapshot/'+ f[:-4] + '.eps')
    
    
print "\a"
