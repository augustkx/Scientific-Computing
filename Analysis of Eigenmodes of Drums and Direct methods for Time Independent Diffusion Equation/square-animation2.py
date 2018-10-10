import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import colors
import numpy as np
from numpy import linalg as LA
import math
from matplotlib import animation

def column(matrix, i):#extract columns from a matrxi
    return [row[i] for row in matrix]


c=1
Lx=1
Ly=Lx=1#square
dx=0.02
dy=dx
Ix=int(Lx/dx+1)#number of points 
Iy=int(Ly/dy+1)#number of points 


#construct the matrix
v = np.zeros((Ix*Iy,Ix*Iy))

for col in range(Ix-2):
    for j in range(Iy*(col+1)+1,Iy*(col+1)+Iy-1):
        v[j,j-Iy]=1/(dx)**2
        v[j,j-1]=1/(dx)**2
        v[j,j]=-4/(dx)**2
        v[j,j+1]=1/(dx)**2
        v[j,j+Iy]=1/(dx)**2

    

eigenvalue, vector = LA.eig(v)        
#eigenvalue, vector = sparse.linalg.eigs(v)   
eigenvalue, vector =np.real(eigenvalue), np.real(vector)


print(eigenvalue)# eigenvalues



K1 = max(i for i in eigenvalue if i < 0)
V1= column(vector, eigenvalue.tolist().index(K1))

K2= max(i for i in eigenvalue if i < 0 and i != K1)
V2= column(vector, eigenvalue.tolist().index(K2))

K3= max(i for i in eigenvalue if i < 0 and i != K1 and i != K2)
V3= column(vector, eigenvalue.tolist().index(K3))

K4= max(i for i in eigenvalue if i < 0 and i != K1 and i != K2 and i != K3 )
V4= column(vector, eigenvalue.tolist().index(K4))

K5= max(i for i in eigenvalue if i < 0 and i != K1 and i != K2 and i != K3 and i != K4 )
V5= column(vector, eigenvalue.tolist().index(K5))
v=V5
v = np.mat(v)
v=v.reshape(Iy,Ix)#matrix


A=1
B=1
T=2000

u = np.zeros((T,Ix,Iy)) #initialize the solution matrix

s=0
for t in [i/1000 for i in range(1,T)]:
     s+=1
     T_separation=A*math.cos(c*(-K1)**0.5*t)+B*math.sin(c*(-K1)**0.5*t)

     for j in range(Iy):
         for i in range(Ix):
              u[s,j,i]=T_separation*v[j,i]
     print(t)

fig = plt.figure()
ax = plt.axes() #xlim=(0, 10), ylim=(0, 10))
     
cmap = mpl.cm.jet
bounds=np.linspace(-0.05, 0.05, num=21)
norm = colors.BoundaryNorm(bounds, cmap.N)
   
im = plt.imshow(u[1,:,:], interpolation='none', origin='lower',cmap=cmap, norm=norm) 
#
##show a color bar
plt.colorbar(im, cmap=cmap, norm=norm, boundaries=bounds, ticks=bounds)

#a text label, used later to show frame number
txt = plt.text(1,1,'')
#
# function for creating an empty plot
def init():
    im.set_data(u[1,:,:])
    return [im]

# i is the frame number
def animate(i,fig,im):
    a=im.get_array()
    a=u[i,:,:]
    im.set_array(a)
    txt.set_text(i)
    return[im, txt]


ani=animation.FuncAnimation(fig,animate,init_func=init,fargs=(fig,im),frames=T)
ani.save('animation-square_f='+str(  '%0.4f' % (  (-K5)**0.5/(2*math.pi),     )      )+'.mp4', fps=20)  

plt.show()

















