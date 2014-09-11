#!/usr/bin/python


import numpy as np
import matplotlib.pyplot as plt
import math
import sys
import re
from decimal import *
import os
import matplotlib as mpl
from scipy.optimize import curve_fit

font = {'family' : 'monospace'}
mpl.rc('font', **font)

pathToGraphFiles ='readyToPlot/'
#graphFiles = [os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser(pathToGraphFiles)) for f in fn  if f.find('.txt') > 0]
graphFiles = [ f for f in os.listdir(pathToGraphFiles) if f.find('.txt') > 0 ]
#print graphFiles

for f in graphFiles:
    opName = f[f.find('OP_')+3:-4]
    pieces = f.split('_')
    graphDescription = ''
    
    for p in range(1, len(pieces)-3, 2):
        tempStr = pieces[p+1]
        tempStr = tempStr.rjust(7 - len(pieces[p])) + '\n'
        graphDescription = graphDescription + pieces[p] + ' = ' + tempStr
    if pieces[len(pieces)-1] != "M.txt":
        pass
        #continue
    #myfile = open(f,'r')
    myfile = open('readyToPlot/'+f,'r')
    xaxis = []
    yaxis = []
    print f
    for line in myfile:
        xaxis.append(float(line[:line.find(' ')]))
        yaxis.append(float(line[line.find(' '):]))
    fig = plt.figure()
#ave = (yaxis[0] + yaxis[1] + yaxis[2])/3.0
 #   if yaxis[0] > 0.5:
#      print yaxis[0]
#continue
    plt.plot(xaxis, yaxis, 'o')
    plt.axes()
 #   plt.xlabel('B')
    plt.xlabel('KbT')
    plt.ylabel(opName)
    plt.minorticks_on
    plt.grid(True, which='both')
    plt.subplots_adjust(right = 0.8)
    if opName[:6] == "Binder":
        plt.ylim(0,0.67)
    if opName[:13] == "Magnetization":
        plt.ylim(0,1)
#popt = []
#pcov = []
#if opName[:6] == "Binder":
#guess = (0.4,0.01,2.0/3.0)#mu, t, n
#plt.plot(xaxis,fermiDirac(np.asarray(xaxis), guess[0], guess[1], guess[2]), 'r-',ls='-')
#popt,pcov = curve_fit(fermiDirac, xaxis, yaxis, guess)
#for p in popt:
#fit = fermiDirac(xaxis, *popt)
#plt.plot(xaxis, fit, 'g-',ls='-')
    plt.figtext(0.81, 0.25, graphDescription)
    #plt.savefig('figures/'+ f[f.find('/'):-4] + '.png')
    plt.savefig('figures/'+ f[:-4] + '.png')
    #plt.savefig('figures/eps/'+ f[:-4] + '.eps')

print "\a"
