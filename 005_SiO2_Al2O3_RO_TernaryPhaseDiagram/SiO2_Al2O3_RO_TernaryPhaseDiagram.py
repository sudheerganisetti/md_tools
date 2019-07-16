#!/usr/bin/python
"""
@Author : S. Ganisetti
@Date   : Fr 17. Aug 22:14:50 CEST 2018

          CaO
          /\
         /  \
        /    \
       /      \
      /________\
    SiO2     Al2O3

m_Si + m_Al + m_Ca = 100

X=(100-m_Si+m_Al)*0.5
Y=sqrt(3/4)*(100-m_Si-m_Al)

"""
import numpy as np
import math
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import subprocess
""" **************** Plot Ternary Phase Diagram **************** """
def plotTernaryPhaseDiagram(x_coord,y_coord,Property1,SampleNames):
  from matplotlib import use
  use('TkAgg')
  min1=min(Property1)
  max1=max(Property1)
  min1=int(min1)
  if int(max1) != max1:
    max1=int(max1)+1
  # a small trick to overcome the singularities if min1=max1
  if min1 == max1:
    min1=min1-0.0001
  diff1=1.0*(max1-min1)
  #Norm_Property1=Property1
  #for i in xrange(len(Property1)):
  #  Norm_Property1[i]=(Property1[i]-min1)/diff1

  fig1= plt.subplot(111)
  plt.xlim(-20,120)
  plt.ylim(-20,100)
  plt.tick_params(axis='both',which='both',top=False,bottom=False,left=False,right=False,labelbottom=False,labelleft=False)
  # triangle outlines
  fig1.plot([0.0,100.0],[0,0],'-k',linestyle='solid')				# Al2O3, bottom axis
  fig1.plot([0.0,50.0],[0,86.6025403784439],'-k',linestyle='solid')		# SiO2, left axis
  fig1.plot([50.0,100.0],[86.6025403784439,0.0],'-k',linestyle='solid')		# CaO, right axis

  #help(mpl.text.Text)
  # Names of the phases
  fig1.text(25.000000-16,43.301270+8,r'$SiO_{2}$',fontsize=20,rotation=60)
  fig1.text(50.000000-9,0.0000000-12,r'$Al_{2}O_{3}$',fontsize=20)
  fig1.text(75.000000+4,43.301270+5,r'$CaO$',fontsize=20,rotation=-60)
  #fig1.text(-7.0,-7.0-1,r'$SiO_{2}$',fontsize=20)
  #fig1.text(91.0,-7.0-1,r'$Al_{2}O_{3}$',fontsize=20)
  #fig1.text(43.0+0.5,90.0,r'$CaO$',fontsize=20)

  # Major Arrows
  #plt.arrow(0.0,-5.0-1,100,0.0,shape='full',head_starts_at_zero=True,head_width=2,length_includes_head=True,color='k',linestyle='solid')                            # bottom
  #plt.arrow(43.0-1,81.8196481400959,-41,-71.0140831103241,shape='full',head_starts_at_zero=True,head_width=2,length_includes_head=True,color='k',linestyle='solid') # left
  #plt.arrow(99.5+1.5,8,-40,69.2820323027552,shape='full',head_starts_at_zero=True,head_width=2,length_includes_head=True,color='k',linestyle='solid')                 # right

  # Guiding lines without arrows
  fig1.plot([37.500000,75.000000],[64.951905,0.000000],'-k',linestyle='dashed') # SiO2, 25 guiding dashed line
  fig1.plot([25.000000,50.000000],[43.301270,0.000000],'-k',linestyle='dashed') # SiO2, 50 guiding dashed line
  fig1.plot([12.500000,25.000000],[21.650635,0.000000],'-k',linestyle='dashed') # SiO2, 75 guiding dashed line
  fig1.plot([25.000000,62.500000],[0.000000,64.951905],'-k',linestyle='dashed') # Al2O3, 25 guiding dashed line
  fig1.plot([50.000000,75.000000],[0.000000,43.301270],'-k',linestyle='dashed') # Al2O3, 50 guiding dashed line
  fig1.plot([75.000000,87.500000],[0.000000,21.650635],'-k',linestyle='dashed') # Al2O3, 75 guiding dashed line
  fig1.plot([12.500000,87.500000],[21.650635,21.650635],'-k',linestyle='dashed') # CaO, 25 guiding dashed line
  fig1.plot([25.000000,75.000000],[43.301270,43.301270],'-k',linestyle='dashed') # CaO, 50 guiding dashed line
  fig1.plot([37.500000,62.500000],[64.951905,64.951905],'-k',linestyle='dashed') # CaO, 75 guiding dashed line

  # Arrows are drawn separately corresponding to each guidline
  plt.arrow(75.000000,0,0.01,-0.015,shape='full',head_starts_at_zero=True,head_width=1,length_includes_head=True,color='k',linestyle='solid')
  plt.arrow(50.000000,0,0.01,-0.015,shape='full',head_starts_at_zero=True,head_width=1,length_includes_head=True,color='k',linestyle='solid')
  plt.arrow(25.000000,0,0.01,-0.015,shape='full',head_starts_at_zero=True,head_width=1,length_includes_head=True,color='k',linestyle='solid')
  plt.arrow(62.500000,64.951905,0.01,0.015,shape='full',head_starts_at_zero=True,head_width=1,length_includes_head=True,color='k',linestyle='solid')
  plt.arrow(75.000000,43.301270,0.01,0.015,shape='full',head_starts_at_zero=True,head_width=1,length_includes_head=True,color='k',linestyle='solid')
  plt.arrow(87.500000,21.650635,0.01,0.015,shape='full',head_starts_at_zero=True,head_width=1,length_includes_head=True,color='k',linestyle='solid')
  plt.arrow(12.500000,21.650635,-0.1,0.0,shape='full',head_starts_at_zero=True,head_width=1,length_includes_head=True,color='k',linestyle='solid')
  plt.arrow(25.000000,43.301270,-0.1,0.0,shape='full',head_starts_at_zero=True,head_width=1,length_includes_head=True,color='k',linestyle='solid')
  plt.arrow(37.500000,64.951905,-0.1,0.0,shape='full',head_starts_at_zero=True,head_width=1,length_includes_head=True,color='k',linestyle='solid')

  # Arrows are drawn separately and outside of triangle corresponding to each guidline
  #plt.arrow(75.000000,0,0.01,-0.015,shape='full',head_starts_at_zero=True,head_width=2,length_includes_head=False,color='k',linestyle='solid')
  #plt.arrow(50.000000,0,0.01,-0.015,shape='full',head_starts_at_zero=True,head_width=2,length_includes_head=False,color='k',linestyle='solid')
  #plt.arrow(25.000000,0,0.01,-0.015,shape='full',head_starts_at_zero=True,head_width=2,length_includes_head=False,color='k',linestyle='solid')
  #plt.arrow(62.500000,64.951905,0.01,0.015,shape='full',head_starts_at_zero=True,head_width=2,length_includes_head=False,color='k',linestyle='solid')
  #plt.arrow(75.000000,43.301270,0.01,0.015,shape='full',head_starts_at_zero=True,head_width=2,length_includes_head=False,color='k',linestyle='solid')
  #plt.arrow(87.500000,21.650635,0.01,0.015,shape='full',head_starts_at_zero=True,head_width=2,length_includes_head=False,color='k',linestyle='solid')
  #plt.arrow(12.500000,21.650635,-0.1,0.0,shape='full',head_starts_at_zero=True,head_width=2,length_includes_head=False,color='k',linestyle='solid')
  #plt.arrow(25.000000,43.301270,-0.1,0.0,shape='full',head_starts_at_zero=True,head_width=2,length_includes_head=False,color='k',linestyle='solid')
  #plt.arrow(37.500000,64.951905,-0.1,0.0,shape='full',head_starts_at_zero=True,head_width=2,length_includes_head=False,color='k',linestyle='solid')

  # tics 
  fig1.text(37.500000-4,64.951905,r'$25$',fontsize=10)
  fig1.text(25.000000-4,43.301270,r'$50$',fontsize=10)
  fig1.text(12.500000-4,21.650635,r'$75$',fontsize=10)
  fig1.text(25.000000,0.0000000-4,r'$25$',fontsize=10)
  fig1.text(50.000000,0.0000000-4,r'$50$',fontsize=10)
  fig1.text(75.000000,0.0000000-4,r'$75$',fontsize=10)
  fig1.text(87.500000+1,21.650635,r'$25$',fontsize=10)
  fig1.text(75.000000+1,43.301270,r'$50$',fontsize=10)
  fig1.text(62.500000+1,64.951905,r'$75$',fontsize=10)

  vmin1=min1
  vmax1=max1
  cmap=plt.cm.jet
  norm=mpl.colors.Normalize(vmin=vmin1,vmax=vmax1)
  ticks1=np.linspace(vmin1,vmax1,5)
  ticks1=np.around(ticks1,decimals=2)
  #fig.scatter(x_coord,y_coord,s=Property1*100,c=Property1,vmin=0,vmax=1,cmap=plt.cm.jet)
  #fig.scatter(x_coord,y_coord,s=100,c=Property1,vmin=0,vmax=1,cmap=plt.cm.jet)
  fig1.scatter(x_coord,y_coord,s=100,c=Property1,cmap=cmap,norm=norm,linewidths=0)
  #fig1.scatter(x_coord,y_coord,s=100,c=Property1,cmap=cmap,vmin=vmin1,vmax=vmax1)
  ax1=plt.axes([0.125,0.1,0.775,0.03]) #left,bottom,width,height
  cb1=mpl.colorbar.ColorbarBase(ax1,cmap=cmap,orientation='horizontal',ticks=ticks1,norm=norm)
  ax1.tick_params(labelsize=20)
  plt.savefig('Fig01_TernaryPhaseDiagramOf_SiO2_Al2O3_CaO.png',bbox_inches='tight',dpi=600)
  for i in xrange(len(x_coord)):
    x1=x_coord[i]
    y1=y_coord[i]
    s1=SampleNames[i]
    fig1.text(x1+2,y1,s1,fontsize=10)
  plt.savefig('Fig02_TernaryPhaseDiagramOf_SiO2_Al2O3_CaO.png',bbox_inches='tight',dpi=600)

  #plt.show()

""" **************** Main Function **************** """  
if __name__=="__main__":
  
  tot_argv=len(sys.argv)
  if tot_argv != 2:
     subprocess.call("sudheer_banner")
     print "************* S. Ganisetti *************"
     print "Error: usage is wrong"
     print "./this_program  SiO2_Al2O3_CaO_PropertyValue.data"
     print "This program is to plot the given property value based on the ternary phase diagram with mol percents of SiO2, Al2O3 and CaO"
     print "The property is plotted from minumum to maximum among the given values and changes from blue to red"
     print "****************************************"
     sys.exit(0)
  INPUTFILE1=str(sys.argv[1])
  data1=np.loadtxt(INPUTFILE1,dtype='str')
  SampleNames=data1.T[0].astype(str)
  SiO2_mol_percent=data1.T[1].astype(float)
  Al2O3_mol_percent=data1.T[2].astype(float)
  Property1=data1.T[4].astype(float)

  x_coord=[0.0 for i in xrange(len(SiO2_mol_percent))]
  y_coord=[0.0 for i in xrange(len(Al2O3_mol_percent))]
  for i in xrange(len(SiO2_mol_percent)):
    #x_coord[i]=Al2O3_mol_percent[i]+(100.0-SiO2_mol_percent[i])/2.0
    #y_coord[i]=1.22474487139159*(100.0-SiO2_mol_percent[i])
    #print "%lf  %lf" %(x_coord[i],y_coord[i])
    x_coord[i]=(100.0-SiO2_mol_percent[i]+Al2O3_mol_percent[i])*0.5
    y_coord[i]=0.866025403784439*(100.0-SiO2_mol_percent[i]-Al2O3_mol_percent[i])
    #print "%lf   %lf" %(x_coord[i],y_coord[i])
  plotTernaryPhaseDiagram(x_coord,y_coord,Property1,SampleNames)
  
