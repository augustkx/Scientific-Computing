
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation



c=1#diffusion constant
dt=0.01#time step
tmin=0#initial time
tmax=1#simulating time
K=int((tmax-tmin)/dt) + 1#number of points on t coordinate

#xmin=0
xmax=1#L=1
I=101#number of points on x coordinates
dx=xmax/(I-1) #intervals  

#inital conditions
def f(x):
    #f=np.sin(2*np.pi*x)
    f=np.sin(5*np.pi*x)
    return f
def f2(x):
    f=np.sin(5*np.pi*x)   
    if 1/5 < x < 2/5:
        return f
    else:
        return 0.0
        
def g(x):
    g=0#The string is at rest at t = 0
    return g    
#boundary coditions
print(K,I)
u = np.zeros((K,I)) #initialize the solution matrix
u[::,0]=0
u[::,I-1]=0

   
#=========================SIMULATION========================================
#initial the shape for nonboundary points
for i in range(1,I-1):    
    u[0,i]=f2(i*dx)#can alter the initial conditions by choosing different expressions of f(x)  or choosing f2(x).
    u[1,i]=u[0,i]+dt*g(i*dx)
#finite difference method applied to nonboundary points
for k in range(1,K-2):
    for i in range(1,I-1):
        u[k+1,i]=2*u[k,i]-u[k-1,i]+c*((dt/dx)**2)*(u[k,i+1]-2*u[k,i]+u[k,i-1])    
    
#===============Plot the results at several times in the same ﬁgure==========

plt.plot([l*dx for l in range(0,I)],u[0], label='t=0')
plt.plot([l*dx for l in range(0,I)],u[5], label='t=0.05')
plt.plot([l*dx for l in range(0,I)],u[10], label='t=0.1')    
plt.plot([l*dx for l in range(0,I)],u[20], label='t=0.2')
plt.plot([l*dx for l in range(0,I)],u[50], label='t=0.5')


plt.legend(loc='upper right')    
plt.xlabel('x')
plt.ylabel('u')    
plt.tight_layout()    
plt.savefig('sin(5pix)_condition.png')
plt.show()

#==========================Animation=========================================
# First set up the figure, the axis, and the plot element we want to animate

fig = plt.figure()    
# set up axes with given limits
ax = plt.axes(xlim=(0, 1), ylim=(-1, 1))
line, = ax.plot([], [], lw=2)  
ax.set_xlabel('x')
ax.set_ylabel('u')
# initialization function: plot what is common to each frame
def init():
    line.set_data([], [])
    return line,

x = [l*dx for l in range(0,I)]

def animate(i):
    print(i)
    y = u[i,:]
    line.set_data(x, y)
    return line,
# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=I, interval=20, blit=True)

# save movie to file. Requires ffmpeg to be installed.
anim.save('animation-sin(5pix)_condition.mp4', fps=20)     
plt.show()


