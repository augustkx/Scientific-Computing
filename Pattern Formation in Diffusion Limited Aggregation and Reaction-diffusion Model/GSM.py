'''Gray-Scott model'''

from matplotlib import colors
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

dt=1#time step
xmin=0
xmax=200#L=1
dx=1#intervals  
ymax=xmax#L=1
dy=dx #intervals  
I=int((xmax-xmin)/dx) + 1 #number of points on x coordinates
d_u=0.16
d_v=0.08
#f,k=0.035,0.06# mu , stripes alter courses to avoid colliding.
#f,k=0.016,0.05# apha, time-dependent and consists of fledgling spirals.
#f,k=0.03,0.057 # thetar,stripes stop when colliding.
f,k=0.0300,0.0630 # eta, time-dependent dots.
#f,k=0.034,0.065#lamda,time-independent, hexagonal grids.


#=============initialize the solution matrix===================================
c_u = np.empty((I+1,I+1)) 
c_u.fill(1)
print(c_u)
c_v = np.zeros((I+1,I+1))
for j in range(97,104):
    for i in range(97,104):
        c_u[j,i]=0.5#inital conditions
        c_v[j,i]=0.25#inital conditions
noise=0#no noise

#==============adding noise====================================================
noise=1
for j in range(0,I-1):
    for i in range(0,I-1):
        if np.random.random()<0.01:
            c_u[j,i] -= 0.5 * np.random.random()#inital conditions
            c_v[j,i] += 0.75 * np.random.random()#inital conditions
            noise=1

#=========================SIMULATION===========================================

#finite difference method applied to nonboundary points
for kk in range(1,(int(1.6E4)+1)):
    c_new_u = np.zeros((I+1,I+1)) 
    c_new_v = np.zeros((I+1,I+1))

    for j in range(1,I):
        for i in range(1,I):#leave the boundary points later to process
           c_new_u[j,i]=c_u[j,i]+( (d_u)*(c_u[j,i+1]+c_u[j,i-1]+c_u[j+1,i]+c_u[j-1,i]-4*c_u[j,i])-c_u[j,i]*c_v[j,i]*c_v[j,i] + f * (1 - c_u[j,i])  )      
           c_new_v[j,i]=c_v[j,i]+( (d_v)*(c_v[j,i+1]+c_v[j,i-1]+c_v[j+1,i]+c_v[j-1,i]-4*c_v[j,i])+c_u[j,i]*c_v[j,i]*c_v[j,i]-(f + k) * c_v[j,i]  )

        c_new_u[j,I]=c_new_u[j,1]#apply the periodic boundary conditions
        c_new_u[j,0]=c_new_u[j,I-1]#apply the periodic boundary conditions
        c_new_v[j,I]=c_new_v[j,1]#apply the periodic boundary conditions
        c_new_v[j,0]=c_new_v[j,I-1]#apply the periodic boundary conditions
    for i in range(0,I+1): 
        c_new_u[0,i]=c_new_u[I-1,i]#apply the periodic boundary conditions
        c_new_u[I,i]=c_new_u[1,i]#apply the periodic boundary conditions 
        c_new_v[0,i]=c_new_v[I-1,i]#apply the periodic boundary conditions
        c_new_v[I,i]=c_new_v[1,i]#apply the periodic boundary conditions 
    
#========================Plotting============================================================================

    if kk in [1,100,int(1E3),int(2E3),int(3E3),int(4E3),int(5E3),int(6E3),int(7E3),int(8E3),int(9E3),int(1E4),int(1.2E4),int(1.4E4),int(1.6E4)]:       
        cmap = mpl.cm.jet
        bounds=np.linspace(0, 1, num=21)
        norm = colors.BoundaryNorm(bounds, cmap.N)
        img = plt.imshow(c_new_u, interpolation='nearest', origin='lower',cmap=cmap, norm=norm)
        plt.gca().invert_yaxis()
        plt.colorbar(img, cmap=cmap, norm=norm, boundaries=bounds, ticks=bounds)
        plt.title('Concentration of U _ t='+str(kk))
        plt.show()
        plt.savefig(str(noise)+'U_f='+str(f)+'_k='+str(k)+'_t=' + str(kk) +'.png')
        plt.clf()  
        cmap = mpl.cm.jet
        bounds=np.linspace(0, 1, num=21)
        norm = colors.BoundaryNorm(bounds, cmap.N)
        img = plt.imshow(c_new_v, interpolation='nearest', origin='lower',cmap=cmap, norm=norm)
        plt.gca().invert_yaxis()
        plt.colorbar(img, cmap=cmap, norm=norm, boundaries=bounds, ticks=bounds)
        plt.title('Concentration of V _ t='+str(kk))
        plt.show()
        plt.savefig(str(noise)+'V_f='+str(f)+'_k='+str(k)+'_t=' + str(kk) +'.png')
        plt.clf()     
   
    #update the matrix
    c_u=c_new_u
    c_v=c_new_v
    

    
    
    
    
    
    
    
    