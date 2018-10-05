import numpy as np
from numpy import log
import matplotlib.pyplot as plt
from collections import defaultdict
import pandas as pd
from sys import setrecursionlimit
setrecursionlimit(15000)


n=256
lambdaa=0.03


#DISTRIBUTIONS
####################################################


#ZiPhz
#s=2
#range(1,n+1)
#F=[i**(-s) for i in range(1,n+1)]
#P=[0]+[i/sum(F) for i in F ]

#EASYTOCHECK
#P=[0]+[0.25,0.25,0.125,0.125,0.0625,0.0625,0.0625,0.0625]

#UNIFORM
#P=[0]+[1.0/n for i in range(n)]
###################################################

#P=[0]+[0.5] +[ 0.5/255 for i in range(255)] 

#P=[0]+[0.4436764149]+[0.0021816611 for i in range(255)] 
P=[0]+[5.0/10]+[ 5.0/(10*255) for i in range(255)]
LE=[0]+[0/n for i in range(n)]
g = np.zeros((n+1,n+1,n+1))
 
  
def ent():
    entr=[-P[i]*log(P[i])/log(2) for i in range(1,n+1)]
    return entr


def best_a_b(a,b):
 g = np.zeros((n+1,n+1,n+1))    
 c = np.zeros((n+1,n+1,n+1))
 LE=[0]+[0/n for i in range(n)]
 for i in range(1,n+1):
    for d in range(1,n+1):
        g[i,i,d]=a*P[i]*d+b*P[i]*d**2
        c[i,i,d]=i

 for j in range(1,n):                
  for i in range(1,n+1-j) :
    for d in [ n-j-1-x for x in range(n-j)]:
         L=[g[i,k-1,d+1]+g[k,i+j,d+1] for k in range(i+1,i+j+1)]
         g[i,i+j,d]=min(L)
         c[i,i+j,d]=range(i+1,i+j+1)[L.index(min(L))]
         #print g[i,i+j,d]

 for i in range(1,n+1):
     LE[i]=0
     x=1
     y=n
     z=int(c[x,y,0]) 
     #print LE, z
     while (i!=x or x!=z-1) and (i!=z or z!=y):     
      if i in range(x,z):
         LE[i]=LE[i]+1
         x=x
         y=z-1
         z=int(c[x,y,LE[i]])
        # print LE, z
      elif i in range (z,y+1):
         LE[i]=LE[i]+1
         x=z
         y=y
         z=int(c[x,y,LE[i]])
         #print LE, z 
     LE[i]=LE[i]+1 
     #print LE
            
                  
         
 return LE         

def EL2(LEN):
    S=[P[i]*LEN[i]**2 for i in range(n+1)] 
    return sum(S)
    
def EL(LEN):
    S=[P[i]*LEN[i] for i in range(n+1)] 
    return sum(S)    

def dele(E,E2, lambdaa):
    return lambdaa*E2/(2*(1-lambdaa*E))+E

def Span(p_1,p_2):
    global pointer,A,B,EL_1,EL_2,D
    a=EL_2[p_1]-EL_2[p_2]
    b=EL_1[p_2]-EL_1[p_1]
    LC=best_a_b(a,b)
    I_1=EL(LC)
    I_2=EL2(LC)
    #print I_1,I_2, pointer,'T'
    A=A+[a]
    B=B+[b]
    EL_1+=[I_1]
    EL_2+=[I_2]
    D+=[dele(I_1,I_2, lambdaa)]
    pointer=pointer+1
    if (a*I_1+b*I_2 == a*EL_1[p_1]+b*EL_2[p_1]) or (a*I_1+b*I_2 == a*EL_1[p_2]+b*EL_2[p_2]):
        print '\n'
        print 'zero' ,(a*I_1+b*I_2), a*EL_1[p_1]+b*EL_2[p_1] ,a*EL_1[p_2]+b*EL_2[p_2]
        return [] 
    else:
      #Span1(E_1,E_2,I_1,I_2)
     # Span1(I_1,I_2,H_1,H_2)
     # print I_1,I_2,pointer,'T'
      return  [(p_1,pointer-1),(pointer-1,p_2)] 
      
     
        
    
                

def main():
 global pointer,A,B,EL_1,EL_2,D  
 global List
 List=[] 
 pointer=0   
 A=[]
 B=[]
 EL_1=[]
 EL_2=[]
 D=[]
 a=1
 b=0
 LC=best_a_b(a,b)
 A=A+[a]
 B=B+[b]
 EL_1+=[EL(LC)]
 EL_2+=[EL2(LC)]
 D+=[dele(EL_1[pointer],EL_2[pointer], lambdaa)]
 pointer=pointer+1
 a=0
 b=1
 LC=best_a_b(a,b)
 A=A+[a]
 B=B+[b]
 EL_1+=[EL(LC)]
 EL_2+=[EL2(LC)]
 D+=[dele(EL_1[pointer],EL_2[pointer], lambdaa)]
 pointer=pointer+1
 #print EL_1,EL_2,pointer
 List=[(0,1)]
 #Temp=Span(List[0][0],List[0][1])
 #List=List[1:]
 #List=List+Temp
 if EL_1[0]==EL_1[1] and EL_2[0]==EL_2[1]:
     return
 else:
    Temp=Span(List[0][0],List[0][1])
    List=List[1:]
    List=List+Temp
    # return   
 while len(List) >=1:
     Temp=Span(List[0][0],List[0][1])
     List=List[1:]
     List=List+Temp
 print 'minimum delay is', min(D) 
 DELAY=[]       
 for lambdaaa in [0.03,0.06,.09,.12,.15]:
     td=[dele(EL_1[i],EL_2[i], lambdaaa) for i in  range(len(EL_2))]
     DELAY+=[min(td)]
 print 'Delay is',DELAY      
main()




