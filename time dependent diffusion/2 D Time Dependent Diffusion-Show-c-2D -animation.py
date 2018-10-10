
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





          
for k in range(1,K):
    c = np.zeros((K,I,I+1)) #initialize the solution matrix
    c[::,0,::]=0#apply the boundary conditions
    c[::,I-1,::]=1#apply the boundary conditions    
    for j in range(1,I-1):
        for i in range(1,I):#leave the boundary points later to process
           c[k,j,i]=c[k-1,j,i]+dt*(d/(dx**2))*(c[k-1,j,i+1]+c[k-1,j,i-1]+c[k-1,j+1,i]+c[k-1,j-1,i]-4*c[k-1,j,i])           
           c[k,I-1,i]=1#applyy the boundary conditions
        c[k,j,0]=c[k,j,I-1]#apply the periodic boudary conditions
        c[k,j,I]=c[k,j,1]#apply the periodic boudary conditions
    print(k)       
        
fig = plt.figure()
ax = plt.axes() #xlim=(0, 10), ylim=(0, 10))
im = plt.imshow(c[1,:,:]) 

#show a color bar
plt.colorbar()

#a text label, used later to show frame number
txt = plt.text(1,1,'')

# function for creating an empty plot
def init():
    im.set_data(c[1,:,:])
    return [im]

# i is the frame number
def animate(i,fig,im):
    a=im.get_array()
    a=c[i,:,:]
    im.set_array(a)
    txt.set_text(i)
    return[im, txt]


ani=animation.FuncAnimation(fig,animate,init_func=init,fargs=(fig,im),frames=K)
ani.save('animation-time-dependent.mp4', fps=20)  
plt.show()

   
