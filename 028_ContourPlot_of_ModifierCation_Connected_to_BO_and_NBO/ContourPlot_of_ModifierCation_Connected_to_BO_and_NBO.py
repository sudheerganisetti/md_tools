import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def scatter_plot(i,j,s,c):
    plt.scatter(i, j, s=s, facecolors=c, edgecolor=c,    label="", linewidth=2)

s0   ='/data2/ganisetti/01_Amar_Sudheer/02_Pedone2006/04_NAP/375N_250A_375P/02_3000atoms/10_properties/5_samples_average_of_BO_and_NBO_connected_to_Na/5_samples_average_of_BO_and_NBO_connected_to_Na.data'
s50  ='/data2/ganisetti/01_Amar_Sudheer/02_Pedone2006/05_NAPS/349N_233A_302P_116S__NAPS050/02_3274_atoms/10_properties/5_samples_average_of_BO_and_NBO_connected_to_Na/5_samples_average_of_BO_and_NBO_connected_to_Na.data'
s100 ='/data2/ganisetti/01_Amar_Sudheer/02_Pedone2006/05_NAPS/326N_218A_239P_217S__NAPS100/01_3074atoms/10_properties/5_samples_average_of_BO_and_NBO_connected_to_Na/5_samples_average_of_BO_and_NBO_connected_to_Na.data'
s150 ='/data2/ganisetti/01_Amar_Sudheer/02_Pedone2006/05_NAPS/306N_204A_184P_306S__NAPS150/02_2902atoms/10_properties/5_samples_average_of_BO_and_NBO_connected_to_Na/5_samples_average_of_BO_and_NBO_connected_to_Na.data'

nap     = np.loadtxt(s0)
naps050 = np.loadtxt(s50)
naps100 = np.loadtxt(s100)
naps150 = np.loadtxt(s150)


plt.rcParams.update({'font.family': 'arial'})
plt.rcParams.update({'font.size': 15})
plt.rcParams.update({'figure.figsize': [5,5]})
fig1=plt.figure()
ax1= fig1.add_subplot(111)

temp1=0
factor1=10.0   # to control the radius
sample_to_color={0:'red',1:'orange',2:'blue',3:'green'}
for i in range(11):
  for j in range(11):
    z1 = nap[temp1][2]
    z2 = naps050[temp1][2]
    z3 = naps100[temp1][2]
    z4 = naps150[temp1][2]
    temp2=np.array([z1,z2,z3,z4])
    temp4=list(temp2)
    temp3=list(temp2)
    temp3.sort()
    for k in [3,2,1]:
      s=temp3[k]*temp3[k]*factor1
      if ( temp3[k] == temp3[k-1] ) and (temp3[k] != 0):
        s=temp3[k]*temp3[k]*(factor1+0.1)
      c=sample_to_color[temp4.index(temp3[k])]
      scatter_plot(i,j,s,c)
    s=temp3[0]*temp3[0]*factor1
    c=sample_to_color[temp4.index(temp3[0])]
    scatter_plot(i,j,s,c)
    temp1=temp1+1


# plotting dummy values to print labels
plt.scatter(0,0,0, facecolors='red', edgecolor="red",    label="Si0",linewidth=1)
plt.scatter(0,0,0, facecolors='orange', edgecolor="orange", label="Si50",linewidth=1)
plt.scatter(0,0,0, facecolors='blue', edgecolor="blue",   label="Si100",linewidth=1)
plt.scatter(0,0,0, facecolors='green', edgecolor="green",  label="Si150",linewidth=1)

ax1.set_xlabel("Number of BO")
ax1.set_ylabel("Number of NBO")
plt.xlim([-0.6, 7.5])
plt.ylim([-0.6, 7.5])
plt.xticks([0,2,4,6])
plt.yticks([0,2,4,6])

legend1 = ax1.legend(markerscale=1,loc=1)
# controlling legend markers size
legend1.legendHandles[0]._sizes = [150]
legend1.legendHandles[1]._sizes = [150]
legend1.legendHandles[2]._sizes = [150]
legend1.legendHandles[3]._sizes = [150]
#fig1.text("NAP",color="red")
plt.savefig('Fig01_Na_in_BO_and_NBO_density_plot.png',bbox_inches='tight',dpi=600)

