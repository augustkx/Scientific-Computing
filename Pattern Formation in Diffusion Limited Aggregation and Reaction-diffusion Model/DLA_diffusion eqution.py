
import matplotlib.pyplot as plt
#from matplotlib import mpl
import numpy as np
from scipy.special import erfc
xmax=1#L=1
I=201#number of points on x coordinates
dx=xmax/(I-1) #intervals  
ymax=1#L=1
dy=ymax/(I-1) #intervals  
eta=1.5#other parameter values: 0,0.5,2,4,10.   the parameters for the probability for growth
w=1#1.8 for eta=0,1,2,4,10 over correction weight.

##======================Analytic solution to set up the empty system========================================================
def analytic(x,t):
    D=1
    N = 30# upper limit for the sum
    c = 0
    d = 2*np.sqrt(D*t)
    for i in range (0,N):
        c += erfc((1-x+2*i) / d) - erfc((1+x+2*i) / d)
    return c
c_y=[analytic(x,1) for x in [l*dy for l in range(0,I)]]# analytic solution
c=np.array([(c_y),]*(I+1)).transpose()#initial concentration domain solved analytically.

aggregation=[]#the aggregatioin of the index belonging to object

aggregation.append([0,int(I/2)])
c[0,int(I/2)]=0                              #a single seed,i.e.sink
iteration=0

while iteration <600:    #the number of chances for the object to grow
    iteration+=1

    candidate=[]
    pro_denominator_each=[]
    for j in range(0,I):
        for i in range(1,I+1):
            if ([j,i] not in aggregation) and ( ([j+1,i] in aggregation) or ([j-1,i] in aggregation) or ([j,i+1] in aggregation) or([j,i-1] in aggregation)):
                candidate.append([j,i])
                pro_denominator_each.append((c[j,i])**eta)
    pro_denominator=np.sum(pro_denominator_each)
    for index in candidate:
        r=np.random.random()
        if r < ((c[index[0],index[1]])**eta)/pro_denominator :#growth probability
            aggregation.append(index)#include the grid into the object
            c[index[0],index[1]]=0#sink
            

    #================================inplement SOR simulation===============================================
    
    tolerance=0.00001# stopping criterion
    sigma=100
    s=0#the number of iterations of SOR method in a single interval between the growth steps.
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
#============delete it when handing up=================================================        
    if  iteration==100 or iteration==200 or iteration==300 or iteration==400 or iteration==500 :   
        np.savetxt('c_eta='+str(eta)+'.txt',c)
        for index in aggregation:
            c[index[0],index[1]]=-0.2#show the object more clearly
        np.savetxt('c_changed_eta='+str(eta)+'.txt',c)    
        for j in range(0,I):
            for i in range(1,I+1):
                if c[j,i]==-0.2:
                    aggregation.append([j,i])
                    c[j,i]=0 
#======================================================================================

     
#np.savetxt('c_eta='+str(eta)+'.txt',c)
for index in aggregation:
    c[index[0],index[1]]=-0.2#show the object more clearly
#np.savetxt('c_changed_eta='+str(eta)+'.txt',c)
plt.imshow(c,cmap='gist_rainbow')  
plt.gca().invert_yaxis()
plt.colorbar()
plt.xlim(0,I)
plt.ylim(0,I)
plt.savefig('DLA_eta='+str(eta)+'.png')  
plt.show()  
    
    