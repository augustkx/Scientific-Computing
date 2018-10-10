
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import erfc
xmax=1#L=1
I=201#number of points on x coordinates
dx=xmax/(I-1) #intervals  
ymax=1#L=1
dy=ymax/(I-1) #intervals  
eta=1#the parameters for the probability for growth

#======================Analytic solution to set up the empty system========================================================
def analytic(x,t):
    D=1
    N = 30                                # upper limit for the sum
    c = 0
    d = 2*np.sqrt(D*t)
    for i in range (0,N):
        c += erfc((1-x+2*i) / d) - erfc((1+x+2*i) / d)
    return c
c_y=[analytic(x,1) for x in [l*dy for l in range(0,I)]]# analytic solution
c=np.array([(c_y),]*(I+1)).transpose()

#=======================simulate the growth model=========================================================================
weight=[i/10 for i in range(12,20)]
ite_con_aver=[]#the average number of iterations SOR method used  for each w parameter.
ite_con_sd=[]
for w in weight:
    ite_convergence=[]#the number of iterations SOR method used  for this w parameter in each PED simulation.
    #a single seed
    c[0,int(I/2)]=0 
    iteration=0
    aggregation=[]#the aggregatioin of the object
    aggregation.append([0,int(I/2)])    
    while iteration <100:    
        iteration+=1

        candidate=[]
        pro_enominator_each=[]
        for j in range(0,I):
            for i in range(1,I+1):
                if ([j,i] not in aggregation) and ( ([j+1,i] in aggregation) or ([j-1,i] in aggregation) or ([j,i+1] in aggregation) or([j,i-1] in aggregation)):
                    candidate.append([j,i])
                    pro_enominator_each.append((c[j,i])**eta)
        pro_denominator=np.sum(pro_enominator_each)
        for index in candidate:
            r=np.random.random()
            if r < ((c[index[0],index[1]])**eta)/pro_denominator :#growth probability
               aggregation.append(index)#include the grid into the object
               c[index[0],index[1]]=0#sink
        
        #================================inplement SOR simulation===============================================
       
        tolerance=0.00001# stopping criterion
        sigma=100
        s=0
        
        while sigma >= tolerance:
            sigma=0
            
        
            for j in range(1,I-1):
                for i in range(1,I):#leave the periodic boundary points later to process
                    if ([j,i] not in aggregation):#if not part of object
                        if abs( w*(c[j,i+1]+c[j,i-1]+c[j+1,i]+c[j-1,i])/4 + (1-w)* c[j,i]-c[j,i] ) > sigma:
                            sigma=abs(w*(c[j,i+1]+c[j,i-1]+c[j+1,i]+c[j-1,i])/4 + (1-w)* c[j,i]-c[j,i] )            
                        c[j,i]=w*(c[j,i+1]+c[j,i-1]+c[j+1,i]+c[j-1,i])/4 + (1-w)* c[j,i]          
        
                c[j,0]=c[j,I-1]#apply the periodic boudary conditions
                c[j,I]=c[j,1]#apply the periodic boudary conditions
        
            s += 1
            print(iteration,s,sigma)
        ite_convergence.append(s)#the number of iterations SOR method used for each PDE simulation   
    ite_con_aver.append(np.average(ite_convergence))#the average number of iterations SOR method used for each w parameter.  
    ite_con_sd.append(np.std(ite_convergence))
#=======================plot the average number of iterations of SOR method during the object growth for each w parameter.
plt.plot(weight,ite_con_aver,'o-')
plt.errorbar(weight,ite_con_aver,yerr = ite_con_sd)
plt.title('Convergence analysis for different w parameters')
plt.xlabel('w')
plt.ylabel('average number of iterations') 
plt.savefig('Convergence_w.png')
plt.show()  



    
    
    
    
    
    
    
    
    
    
    