
import matplotlib.pyplot as plt
import numpy as np
#from matplotlib import animation


xmax=1#L=1
I=61#number of points on x coordinates
dx=xmax/(I-1) #intervals  
ymax=1#L=1
dy=ymax/(I-1) #intervals  




# initialize the array/matrix for plotting

weight=[i/10 for i in range(10,20)]
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
            print(w,iteration_SOR,sigma)
        iteration_SOR=np.log10(iteration_SOR)
        matrix_iteration_SOR[q,tol-3]=iteration_SOR
        q+=1
    



for q in range(len(weight)):
    plt.plot(range(3,9), matrix_iteration_SOR[q],'o-',label="w="+str(weight[q]))
plt.legend()
plt.title('Convergence analysis for SOR_N='+ str(I-1))
plt.xlim(2,10)
plt.xlabel('p')
plt.ylabel('iteration') 
plt.savefig('Convergence_optimal_N_'+ str(I-1)+'.png')
plt.show()  




