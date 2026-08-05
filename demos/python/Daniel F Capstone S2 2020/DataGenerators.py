import random
from collections import Counter
import numpy as np 

starts=[0, 0.1,0.2] 
ends=[1,1.2,1.3]

m=[0.4,0.5,0.8] #Means
s=[0.2, 0.3, 0.25]#Standard deviations

def GenerateSquareWaveData(ab_norm, start_points=starts,end_points=ends, size=1000):
    x=[] #these are not the same as the x-axis and y-axis that will be plotted later
    y=[]     
    a=0
    list_of_sizes=[0]*len(starts)
    many_samples = random.choices(list(range(len(ab_norm))), weights=ab_norm, k=size)
    dict1=Counter(many_samples)
    for key, item in dict1.items():
        list_of_sizes[int(key)]= item
    while a<len(starts):
        i=0
        while i<list_of_sizes[a]: #ensuring our cts and dsc datasets have same len
            x.append(a)
            y.append(random.uniform(start_points[a],end_points[a]))
            i=i+1
        a=a+1    
    return x,y, len(list_of_sizes)  

#1. Generating continuous data_____________________________
def GenerateSingleGaussianData (ab_norm, list_of_mu =m, list_of_sigma=s,size=1000):    
    x=[] #these are not the same as the x-axis and y-axis that will be plotted later
    y=[] 
    a=0    
    list_of_sizes=[0]*len(list_of_mu)
    many_samples = random.choices(list(range(len(ab_norm))), weights=ab_norm, k=size)
    dict1=Counter(many_samples)
    for key, item in dict1.items():
        list_of_sizes[int(key)]= item
    while a<len(list_of_mu):
        i=0            
        while i<list_of_sizes[a]: 
            x.append(a)
            y.append(random.gauss(list_of_mu[a], list_of_sigma[a]))
            i=i+1     
        a=a+1         

    return x,y, len(list_of_sizes)  


def GenerateUncorrelatedGaussianData(ab_norm, list_of_mu =m, list_of_sigma=s,size=1000):    
    x=[] #these are not the same as the x-axis and y-axis that will be plotted later
    y=[] 
    z=[]
    a=0    
    list_of_sizes=[0]*len(list_of_mu)
    many_samples = random.choices(list(range(len(ab_norm))), weights=ab_norm, k=size)
    dict1=Counter(many_samples)
    for key, item in dict1.items():
        list_of_sizes[int(key)]= item
    while a<len(list_of_mu):
        i=0            
        while i<list_of_sizes[a]: 
            x.append(a)
            y.append(random.gauss(list_of_mu[a], list_of_sigma[a]))
            z.append(random.gauss(0,1))
            i=i+1     
        a=a+1         

    return x,y,z, len(list_of_sizes)

def GenerateCorrelatedGaussianData (ab_norm,size,list_of_mu =m, list_of_sigma=s):    
    x=[] #these are not the same as the x-axis and y-axis that will be plotted later
    y=[]   
    z=[]
    a=0    
    list_of_sizes=[0]*len(list_of_mu)
    many_samples = random.choices(list(range(len(ab_norm))), weights=ab_norm, k=size)
    dict1=Counter(many_samples)
    for key, item in dict1.items():
        list_of_sizes[int(key)]= item
    while a<len(list_of_mu):
        i=0            
        while i<list_of_sizes[a]: #ensuring our cts and dsc datasets have same len
            x.append(a)
            mean=[0, list_of_mu[a]]
            cov=[[1,np.sqrt(list_of_sigma[a])],[np.sqrt(list_of_sigma[a]),list_of_sigma[a]]]
            obs=np.random.multivariate_normal(mean,cov)
            y.append(obs[1])
            z.append(obs[0])
            i=i+1  
        a=a+1       
    
    return x,y,z,len(list_of_sizes)

