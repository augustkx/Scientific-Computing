import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg as LA

def column(matrix, i):#extract columns from a matrxi
    return [row[i] for row in matrix]
def y(index):#return to the y-coordinate position in original domain for M matrix.
    return index % Iy
def x(index):#return to the x-coordinate position in original domain for M matrix.
    return index // Iy



dx=dy=0.05
Ly=Lx=4
Ix=int(Lx/dx+1)#number of points 
Iy=int(Ly/dy+1)#number of points 
print(Ix*Iy)
#============================construct the matrix M===========================================================
m = np.zeros((Ix*Iy,Ix*Iy))


for j in range(Iy+1,Ix*Iy-Iy-1):         
        m[j,j-Iy]=1/(dx)**2
        m[j,j-1]=1/(dx)**2
        m[j,j]=-4/(dx)**2
        m[j,j+1]=1/(dx)**2
        m[j,j+Iy]=1/(dx)**2

#update the rows to zero, whose corresponding points do not belong to the circle domian.     
for i in range(Ix*Iy):
    if  ((x(i)-(Ix-1)/2)**2+(y(i)-(Iy-1)/2)**2)**0.5 >= (Ix-1)/2:# outside the domain. Center:(  (Ix-1)/2,(Iy-1)/2  )
        m[i]=np.zeros(Ix*Iy)
        m[i,i]=1
#update the rows corresponding to the sourse
m[int(2.6/dx)* Iy +3.2/dy]=np.zeros(Ix*Iy)
m[int((2.6/dx)* Iy +3.2/dy),int((2.6/dx)* Iy +3.2/dy)]=1
print((2.6/dx)* Iy +3.2/dy)
#============================construct b vector================================================================
b=np.zeros(Ix*Iy)
b[int((2.6/dx)* Iy +3.2/dy)]=1
#============================solve the algebraic equations=====================================================
u = LA.solve(m, b)

#========================Plotting================================================================================
u = np.mat(u)
u=u.reshape(Iy,Ix)#matrix
plt.imshow(u)
plt.colorbar()
plt.gca().invert_yaxis()
plt.savefig('Direct methods for solving steady state problems.png')

    













