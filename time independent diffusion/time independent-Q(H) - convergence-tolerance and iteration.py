
import matplotlib.pyplot as plt
import numpy as np
#from matplotlib import animation


xmax=1#L=1
I=51#number of points on x coordinates
dx=xmax/(I-1) #intervals  
ymax=1#L=1
dy=ymax/(I-1) #intervals  




# initialize the array/matrix for plotting
array_iteration_GS=[]
array_iteration_J=[]
weight= [0.5,1.5,1.8]#SOR method
matrix_iteration_SOR=np.zeros((len(weight),len(range(3,9))))

for tol in range(3,9):# simulate for various tolerance level
    tolerance=10**(-tol)# stopping criterion
    print(tolerance,"-----------------------")
    
    #=========================SIMULATION-SOR Iteration========================================

    
    q=0
    for w in weight:    
        c = np.zeros((I,I+1)) #initialize the solution matrix
        c[0,::]=0#boundary coditions
        c[I-1,::]=1#boundary coditions       
#        w=1.8#over correction weight
        sigma=100
        iteration_SOR=0
        while sigma >= tolerance:
            sigma=0
          
            for j in range(1,I-1):
                for i in range(1,I):#leave the boundary points later to process
                    if abs( w*(c[j,i+1]+c[j,i-1]+c[j+1,i]+c[j-1,i])/4 + (1-w)* c[j,i]-c[j,i] ) > sigma:
                        sigma=abs(w*(c[j,i+1]+c[j,i-1]+c[j+1,i]+c[j-1,i])/4 + (1-w)* c[j,i]-c[j,i] )            
                    c[j,i]=w*(c[j,i+1]+c[j,i-1]+c[j+1,i]+c[j-1,i])/4 + (1-w)* c[j,i]          
        
                c[j,0]=c[j,I-1]#apply the periodic boudary conditions
                c[j,I]=c[j,1]#apply the periodic boudary conditions
        
            iteration_SOR += 1
            print("SOR",iteration_SOR,sigma)
        iteration_SOR=np.log10(iteration_SOR)
        matrix_iteration_SOR[q,tol-3]=iteration_SOR
        q+=1
    
    #=========================SIMULATION-Gauss-Seidel Iteration========================================
    c = np.zeros((I,I+1)) #initialize the solution matrix
    c[0,::]=0#boundary coditions
    c[I-1,::]=1#boundary coditions
    sigma=100
    iteration_GS=0
    while sigma >= tolerance:
        sigma=0
        for j in range(1,I-1):
            for i in range(1,I):#leave the boundary points later to process
                if abs((c[j,i+1]+c[j,i-1]+c[j+1,i]+c[j-1,i])/4-c[j,i]) > sigma:
                    sigma=abs((c[j,i+1]+c[j,i-1]+c[j+1,i]+c[j-1,i])/4-c[j,i])            
                c[j,i]=(c[j,i+1]+c[j,i-1]+c[j+1,i]+c[j-1,i])/4          
    
            c[j,0]=c[j,I-1]#apply the periodic boudary conditions
            c[j,I]=c[j,1]#apply the periodic boudary conditions
    
        iteration_GS += 1
        print("G",iteration_GS,sigma)
    iteration_GS=np.log10(iteration_GS)
    array_iteration_GS.append(iteration_GS)
       
    #=========================SIMULATION-Jacobi Iteration========================================
    c = np.zeros((I,I+1)) #initialize the solution matrix
    c[0,::]=0#boundary coditions
    c[I-1,::]=1#boundary coditions
    sigma=100
    iteration_J=0
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
        iteration_J += 1
        print("J",iteration_J,sigma)
    iteration_J=np.log10(iteration_J)
    array_iteration_J.append(iteration_J)

plt.plot(range(3,9),array_iteration_J,'o-',label="Jacobi")
plt.plot(range(3,9),array_iteration_GS,'o-',label="JGauss-Seidel")
for q in range(len(weight)):
    plt.plot(range(3,9), matrix_iteration_SOR[q],'^-',label="SOR-w="+str(weight[q]))
plt.legend()
plt.title('Convergence analysis for tolerance level and iteration number')
plt.xlim(2,10)
plt.xlabel('p')
plt.ylabel('iteration') 
plt.savefig('Convergence_3method.png')
plt.show()  




