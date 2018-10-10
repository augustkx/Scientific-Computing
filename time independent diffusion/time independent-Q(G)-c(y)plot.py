
import matplotlib.pyplot as plt
import numpy as np
#from matplotlib import animation


xmax=1#L=1
I=51#number of points on x coordinates
dx=xmax/(I-1) #intervals  
ymax=1#L=1
dy=ymax/(I-1) #intervals  


#=========================SIMULATION-SOR Iteration========================================
c = np.zeros((I,I+1)) #initialize the solution matrix
c[0,::]=0#boundary coditions
c[I-1,::]=1#boundary coditions
w=1.8#over correction weight
tolerance=0.00001# stopping criterion
sigma=100
s=0
while sigma >= tolerance:
    sigma=0
    
    for j in range(1,I-1):
        for i in range(1,I):#leave the boundary points later to process
            if abs( w*(c[j,i+1]+c[j,i-1]+c[j+1,i]+c[j-1,i])/4 + (1-w)* c[j,i]-c[j,i] ) > sigma:
                sigma=abs(w*(c[j,i+1]+c[j,i-1]+c[j+1,i]+c[j-1,i])/4 + (1-w)* c[j,i]-c[j,i] )            
            c[j,i]=w*(c[j,i+1]+c[j,i-1]+c[j+1,i]+c[j-1,i])/4 + (1-w)* c[j,i]          

        c[j,0]=c[j,I-1]#apply the periodic boudary conditions
        c[j,I]=c[j,1]#apply the periodic boudary conditions

    s += 1
    print(s,sigma)
#==============Plot of c(y)===============================================================
c_y=c[::,1]
plt.plot([l*dy for l in range(0,I)],c_y,label="w="+ str(w))
plt.title('SOR Iteration Method')
plt.xlabel('y')
plt.ylabel('c') 
plt.savefig('c(y)_SOR.png')
plt.show()  
 
##=============Exploration of the numerical error of the Gauss-Seidel Iteration Method: the dependence of c value on x.
for x in range(1,I+1):
    plt.plot([l*dy for l in range(0,I)],c[::,x]-c[::,x-1])

plt.title('SOR Iteration Method Error')
plt.xlabel('y')
plt.ylabel('c') 
plt.savefig('c(y)_SOR_error.png')
plt.show()  



#=========================SIMULATION-Gauss-Seidel Iteration========================================
c = np.zeros((I,I+1)) #initialize the solution matrix
c[0,::]=0#boundary coditions
c[I-1,::]=1#boundary coditions

tolerance=0.00001# stopping criterion
sigma=100
s=0
while sigma >= tolerance:
    sigma=0
    for j in range(1,I-1):
        for i in range(1,I):#leave the boundary points later to process
            if abs((c[j,i+1]+c[j,i-1]+c[j+1,i]+c[j-1,i])/4-c[j,i]) > sigma:
                sigma=abs((c[j,i+1]+c[j,i-1]+c[j+1,i]+c[j-1,i])/4-c[j,i])            
            c[j,i]=(c[j,i+1]+c[j,i-1]+c[j+1,i]+c[j-1,i])/4          

        c[j,0]=c[j,I-1]#apply the periodic boudary conditions
        c[j,I]=c[j,1]#apply the periodic boudary conditions

    s += 1
    print(s,sigma)
#==============Plot of c(y)===============================================================
c_y=c[::,1]
plt.plot([l*dy for l in range(0,I)],c_y)
plt.title('Gauss-Seidel Iteration Method')
plt.xlabel('y')
plt.ylabel('c') 
plt.savefig('c(y)_Gauss-Seidel.png')
plt.show()
plt.clf()  
#=============Exploration of the numerical error of the Gauss-Seidel Iteration Method: the dependence of c value on x.

for x in range(1,I+1):
    plt.plot([l*dy for l in range(0,I)],c[::,x]-c[::,x-1])

plt.title('Gauss-Seidel Iteration Method Error')
plt.xlabel('y')
plt.ylabel('c') 
plt.savefig('c(y)_Gauss-Seidel_error.png')
plt.show()  


#=========================SIMULATION-Jacobi Iteration========================================
c = np.zeros((I,I+1)) #initialize the solution matrix
c[0,::]=0#boundary coditions
c[I-1,::]=1#boundary coditions

tolerance=0.00001# stopping criterion
sigma=100
s=0
while sigma >= tolerance:
    sigma=0
    c_new = np.zeros((I,I+1)) 
    c_new[0,::]=0#apply the boundary conditions
    c_new[I-1,::]=1  #apply the boundary conditions       
    for j in range(1,I-1):
        for i in range(1,I):#leave the boundary points later to process
            c_new[j,i]=(c[j,i+1]+c[j,i-1]+c[j+1,i]+c[j-1,i])/4          
          
           # sigma=abs(c_new[j,i]-c[j,i])
            if abs(c_new[j,i]-c[j,i]) > sigma:
                sigma=abs(c_new[j,i]-c[j,i])

        c_new[j,0]=c_new[j,I-1]#apply the periodic boudary conditions
        c_new[j,I]=c_new[j,1]#apply the periodic boudary conditions
    
    c=c_new#update the 
    s += 1
    print(s,sigma)
#print(sigma)
c_y=c[::,1]
plt.plot([l*dy for l in range(0,I)],c_y)
plt.title('Jacobi Iteration Method')
plt.xlabel('y')
plt.ylabel('c') 
plt.savefig('c(y)_Jacobi.png')
plt.show()  




