import matplotlib.pyplot as plt
#from matplotlib import mpl
import numpy as np
from numpy import linalg as LA
import math
c=1
Lx=1

#dx=0.02
#dy=dx
Ix=51#number of points 
Iy=51#number of points 


def column(matrix, i):#extract columns from a matrxi
    return [row[i] for row in matrix]



s=0
for dx in [0.02,0.04,0.0625,0.1]:#spectrum of frequencies with variious number of discretization steps.
#for Lx in  [1,1.5,2]:# spectrum of eigenfrequencies with different size L.
    Ly=Lx#square
    s+=1
    dy=dx
    Ix=int(Lx/dx+1)#number of points 
    Iy=int(Ly/dy+1)#number of points 
#    dx=Lx/(Ix-1)
#    dy=dx
    
    
    #construct the matrix
    v = np.zeros((Ix*Iy,Ix*Iy))
    
    for col in range(Ix-2):
        for j in range(Iy*(col+1)+1,Iy*(col+1)+Iy-1):
        #for j in range(Iy+1,Ix*Iy-Iy-1):    
            
            v[j,j-Iy]=1/(dx)**2
            v[j,j-1]=1/(dx)**2
            v[j,j]=-4/(dx)**2
            v[j,j+1]=1/(dx)**2
            v[j,j+Iy]=1/(dx)**2
    
        
    
    eigenvalue, vector = LA.eig(v)        
    #eigenvalue, vector = sparse.linalg.eigs(v)   
#    eigenvalue, vector =np.real(eigenvalue), np.real(vector)
    
    #print(Ix*Iy)
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
    #==================================Plot the eigenvectors v for some of the smallest eigenvalues=========================
    plt.plot([l for l in range(0,Ix*Iy)], V1, label= 'f='+str('%0.4f' % (       (-K1)**0.5/(2*math.pi),     )))
    plt.plot([l for l in range(0,Ix*Iy)], V2, label= 'f='+str('%0.4f' % (       (-K2)**0.5/(2*math.pi),     )))
    plt.plot([l for l in range(0,Ix*Iy)], V3, label= 'f='+str('%0.4f' % (       (-K3)**0.5/(2*math.pi),     )))
    plt.plot([l for l in range(0,Ix*Iy)], V4, label= 'f='+str('%0.4f' % (       (-K4)**0.5/(2*math.pi),     )))
    plt.plot([l for l in range(0,Ix*Iy)], V5, label= 'f='+str('%0.4f' % (       (-K5)**0.5/(2*math.pi),     )))
    plt.legend(loc='upper right')  
    plt.savefig('vector_square2.png')
    plt.clf()
    
    #================================= spectrum of eigenfrequencies depend on the size L=1 (comment previous plots code before this plotting)==================
#    color=['','ro','bo','go','yo']
#    plt.plot([(-i)**0.5/(2*math.pi) for i in eigenvalue],color[s],label='L='+str(Ly))
#    plt.plot([(-i)**0.5/(2*math.pi) for i in eigenvalue],color[s],label='dx='+str(dx))
#    plt.plot([(-i)**0.5/(2*math.pi) for i in eigenvalue],color[s],label='L='+str(Lx)+', dx='+str(dx))
#plt.title("the spectrum of eigenfrequencies_Square")
#plt.legend(loc='upper right')  
#plt.savefig('eigenfrequencies_square3.png')












