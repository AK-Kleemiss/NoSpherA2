from matplotlib import pyplot as plt
from matplotlib import ticker
import os
import math
import matplotlib.font_manager as font_manager

xray = False
ED = True

font = font_manager.FontProperties(style='normal', size=9)

def add_subplot_axes(ax,rect):
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)    
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]
    subax = fig.add_axes([x,y,width,height])
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax

if xray == True:
  thakkars = [[] for x in range(6)]
  k_vectors = []
  file = open("sfacs.dat","r")
  for line in file.readlines():
      values = line.split(" ")
      k_vectors.append(float(values[0]))
      for i in range(6):
          thakkars[i].append(float(values[i+1]))
  
  print("Reading done!")
  fig = plt.figure()
  axes = fig.add_subplot(1,1,1)
  axes.set_xlabel(r'$k /\AA^{-1}$')
  axes.set_ylabel(r'$|f| /e$')
  if i > 0:
      axes.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
  axes.set_xlim(0,2.0)
  axes.scatter(k_vectors,thakkars[3],s=0.1,linewidths=0.9,facecolors='r',edgecolors='r',label="O (epoxide)")
  axes.plot(k_vectors,thakkars[0],label="O")
  axes.plot(k_vectors,thakkars[1],label="O-")
  axes.plot(k_vectors,thakkars[2],label="O+")
  
  subax = add_subplot_axes(axes,[0.3,0.45,0.68,0.4])
  subax.set_xlim(0.5,1.0)
  subax.set_ylim(1.0,3.0)
  subax.plot(k_vectors,thakkars[0])
  subax.plot(k_vectors,thakkars[1])
  subax.plot(k_vectors,thakkars[2])
  subax.scatter(k_vectors,thakkars[3],s=0.1,linewidths=0.9,facecolors='r',edgecolors='r')
  axes.legend(ncol=4,prop=font)
  fig.savefig("sfacs.png",dpi=600,transparent=True,bbox_inches='tight')
  
  fig = plt.figure()
  axes = fig.add_subplot(1,1,1)
  axes.set_xlabel(r'$k /\AA^{-1}$')
  axes.set_ylabel(r'$f /e$')
  if i > 0:
      axes.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
  axes.set_xlim(0,2.0)
  axes.plot(k_vectors,thakkars[0],label="O")
  axes.plot(k_vectors,thakkars[1],label="O-")
  axes.plot(k_vectors,thakkars[2],label="O+")
  
  subax = add_subplot_axes(axes,[0.3,0.45,0.68,0.4])
  subax.set_xlim(0,0.5)
  subax.plot(k_vectors,thakkars[0])
  subax.plot(k_vectors,thakkars[1])
  subax.plot(k_vectors,thakkars[2])
  
  axes.legend(ncol=3,prop=font)
  fig.savefig("sfacs_spherical_only.png",dpi=600,transparent=True,bbox_inches='tight')
  
  fig = plt.figure()
  axes = fig.add_subplot(1,1,1)
  axes.set_xlabel(r'$k /\AA^{-1}$')
  axes.set_ylabel(r'$f /e$')
  if i > 0:
      axes.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
  axes.set_xlim(0,2.0)
  axes.plot(k_vectors,thakkars[0],label="O")
  axes.plot(k_vectors,thakkars[4],label="O core")
  axes.plot(k_vectors,thakkars[5],label="O valence")
  
  subax = add_subplot_axes(axes,[0.3,0.45,0.68,0.4])
  subax.set_xlim(0,0.5)
  subax.plot(k_vectors,thakkars[0])
  subax.plot(k_vectors,thakkars[4])
  subax.plot(k_vectors,thakkars[5])
  
  axes.legend(ncol=3,prop=font)
  fig.savefig("sfacs_core_val_split.png",dpi=600,transparent=True,bbox_inches='tight')
  
  
  print("finished with X-ray, now ED!")
if ED == True:
  thakkars = [[] for x in range(6)]
  k_vectors = []
  file = open("sfacs_ED.dat","r")
  for line in file.readlines():
      values = line.split(" ")
      k_vectors.append(float(values[0]))
      for i in range(6):
          thakkars[i].append(float(values[i+1]))
  
  print("Reading done!")
  fig = plt.figure()
  axes = fig.add_subplot(1,1,1)
  axes.set_xlabel(r'$k /\AA^{-1}$')
  axes.set_ylabel(r'$|f| /\AA$')
  if i > 0:
      axes.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
  axes.set_xlim(0,2.0)
  axes.set_ylim(-1,3.0)
  axes.scatter(k_vectors,thakkars[3],s=0.1,linewidths=0.9,facecolors='r',edgecolors='r',label="O (epoxide)")
  axes.plot(k_vectors,thakkars[0],label="O")
  axes.plot(k_vectors,thakkars[1],label="O-")
  axes.plot(k_vectors,thakkars[2],label="O+")
  
  subax = add_subplot_axes(axes,[0.3,0.45,0.68,0.4])
  subax.set_xlim(0,0.3)
  subax.set_ylim(-2,6.0)
  subax.plot(k_vectors,thakkars[0])
  subax.plot(k_vectors,thakkars[1])
  subax.plot(k_vectors,thakkars[2])
  subax.scatter(k_vectors,thakkars[3],s=0.1,linewidths=0.9,facecolors='r',edgecolors='r')
  axes.legend(ncol=4,prop=font)
  fig.savefig("sfacs_ED.png",dpi=600,transparent=True,bbox_inches='tight')
  
  fig = plt.figure()
  axes = fig.add_subplot(1,1,1)
  axes.set_xlabel(r'$k /\AA^{-1}$')
  axes.set_ylabel(r'$f /\AA$')
  if i > 0:
      axes.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
  axes.set_xlim(0,2.0)
  axes.set_ylim(-1,3.0)
  axes.plot(k_vectors,thakkars[0],label="O")
  axes.plot(k_vectors,thakkars[1],label="O-")
  axes.plot(k_vectors,thakkars[2],label="O+")
  
  subax = add_subplot_axes(axes,[0.3,0.45,0.68,0.4])
  subax.set_xlim(0,0.5)
  subax.set_ylim(-1,3.0)
  subax.plot(k_vectors,thakkars[0])
  subax.plot(k_vectors,thakkars[1])
  subax.plot(k_vectors,thakkars[2])
  
  axes.legend(ncol=3,prop=font)
  fig.savefig("sfacs_spherical_only_ED.png",dpi=600,transparent=True,bbox_inches='tight')
  
  fig = plt.figure()
  axes = fig.add_subplot(1,1,1)
  axes.set_xlabel(r'$k /\AA^{-1}$')
  axes.set_ylabel(r'$f /\AA$')
  if i > 0:
      axes.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
  axes.set_xlim(0,2.0)
  axes.set_ylim(-1,3.0)
  axes.plot(k_vectors,thakkars[0],label="O")
  axes.plot(k_vectors,thakkars[4],label="O core")
  axes.plot(k_vectors,thakkars[5],label="O valence")
  
  subax = add_subplot_axes(axes,[0.3,0.45,0.68,0.4])
  subax.set_xlim(0,0.5)
  subax.set_ylim(-1,3.0)
  subax.plot(k_vectors,thakkars[0])
  subax.plot(k_vectors,thakkars[4])
  subax.plot(k_vectors,thakkars[5])
  
  axes.legend(ncol=3,prop=font)
  fig.savefig("sfacs_core_val_split_ED.png",dpi=600,transparent=True,bbox_inches='tight')
 
