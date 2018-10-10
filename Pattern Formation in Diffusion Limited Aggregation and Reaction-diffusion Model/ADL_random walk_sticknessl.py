import matplotlib.pyplot as plt
import numpy as np
import random
       
def move(c,direction):
    global position 
    j=position[0]
    i=position[1]
    if direction%2 == 0:#horizontal direction
        if [j,i+direction-3]  not in aggregation:#if walker is not stepping into the clustet,then check whether it moves out of boundary; else, position unchanged.  
            if i+direction-3!= I and i+direction-3!= -1:#not stepping out of boundary
                position=[j,i+direction-3]
            else:#stepping out of boundary
                if i+direction-3== I :
                    position=[j,0]
                else:
                    position=[j,I-1]
    
    else:#vertical direction
        if [j+direction-2,i] not in aggregation:#if walker is not stepping into the clustet,then check whether it moves out of boundary; else, position unchanged. 
            if j+direction-2 != I and j+direction-2 != -1:#if not stepping out of boundary, update the position.Otherwise, unchanged.
                position=[j+direction-2,i]

    
I=300#number of points on x coordinates
c = np.ones((I,I)) 
c[0,int(I/2)]=0 #initial object
stickness=0.75#sticking probability

aggregation=[]#the aggregatioin of the index belonging to object
aggregation.append([0,int(I/2)])


            
size=0
while size < 1000:
    position=[I-1,random.randint(0,I-1)]
    stick=0#0:the walker is free;1:the walker is stuck to the object.
    neighbors=[]
    size+=1#the size of the object.
    for j in range(0,I):
        for i in range(0,I):
            if ([j,i] not in aggregation) and ( ([j+1,i] in aggregation) or ([j-1,i] in aggregation) or ([j,i+1] in aggregation) or([j,i-1] in aggregation)):
                 neighbors.append([j,i])
    while stick == 0:#while the walker still not sticks to the object.
        direction=random.randint(1,4)#1:down;2:left;3:up;4:right.
        #if the direction chosen makes the walker step into the cluster, then produce a direction randomly again until it makes way into an available site.
        while (direction%2 == 0 and ([j,i+direction-3]  in aggregation)) or (direction%2 == 1 and ([j+direction-2,i] in aggregation)):
            direction=random.randint(1,4)
        move(c,direction)
        if position in neighbors:#reach the neighbors of the object.
            r=np.random.random()
            if  r < stickness:
                stick=1
                aggregation.append(position)


#set different values for the cluster points in order to show the cluster image.       
for j in range(0,I):
    for i in range(1,I):
        if [j,i] in aggregation:
            c[j,i]=0        
 
     
#np.savetxt("DLA_random_"+str(stickness)+".txt",c)
plt.imshow(c,cmap='Greys')
plt.gca().invert_yaxis()
plt.colorbar()
plt.xlim(0,I-1)
plt.ylim(0,I-1)
plt.savefig("DLA_random_"+str(stickness)+".png")  
plt.show()  