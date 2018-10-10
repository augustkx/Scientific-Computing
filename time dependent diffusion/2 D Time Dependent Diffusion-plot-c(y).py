
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation


d=1#diffusion constant
dt=0.000025#time step
tmin=0#initial time
tmax=1#simulating time
K=int((tmax-tmin)/dt) + 1#number of points on t coordinate


xmax=1#L=1
I=101#number of points on x coordinates
dx=xmax/(I-1) #intervals  
ymax=1#L=1
dy=ymax/(I-1) #intervals  

#check the stability of the scheme
if 4*dt*d>dx**2:
    print("Pleas reset the parameters!")
    quit()


#boundary coditions & inital conditions
c = np.zeros((I,I+1)) #initialize the storage of solution matrix
c[0,::]=0
c[I-1,::]=1

#=========================SIMULATION========================================

#finite difference method applied to nonboundary points
for k in range(1,K):
    c_y=[0]#initialize the c(y) array.
    c_new = np.zeros((I,I+1)) #initialize the solution matrix
    c_new[0,::]=0#apply the boundary conditions
    c_new[I-1,::]=1  #apply the boundary conditions   
    for j in range(1,I-1):
        for i in range(1,I):#leave the boundary points later to process
           c_new[j,i]=c[j,i]+dt*(d/(dx**2))*(c[j,i+1]+c[j,i-1]+c[j+1,i]+c[j-1,i]-4*c[j,i])           
        c_new[j,0]=c_new[j,I-1]#apply the periodic boudary conditions
        c_new[j,I]=c_new[j,1]#apply the periodic boudary conditions
        
        #plot c(y) for different times:t=0.001,0.01,0.1, and 1.In order to compare with the analytic solutions.
        print(k,j)
        if k == 0.001/dt or k == 0.01/dt  or k == 0.1/dt or k == 1/dt:
            c_y.append(c_new[j,1])

    if k == 0.001/dt or k == 0.01/dt  or k == 0.1/dt or k == 1/dt:
        c_y.append(1)#boundary conditions
        plt.plot([l*dy for l in range(0,I)],c_y,label='t='+str(k*dt))
        plt.legend(loc='upper right')  
        plt.show()  
    
    c=c_new #update the matrix

plt.legend(loc='upper right')          
plt.xlabel('y')
plt.ylabel('c')    
plt.tight_layout()    
plt.savefig('c(y).png')
plt.show()        

 

