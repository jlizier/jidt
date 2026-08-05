from scipy.integrate import quad
import numpy as np

starts=[0, 0.1,0.2] 
ends=[1,1.2,1.3]
list_abund=[2,10,5] # determines how many discrete values & their relative prevalence
weights=list_abund/sum(np.array(list_abund)) #normalizing
    
def mu_multiple(y,weights=weights,starts=starts,ends=ends):
    j=0    
    subtotal=0
    while j<len(weights):        
        length=ends[j]-starts[j]
        ind_result=weights[j]/length         
        if y>starts[j] and y<ends[j]:                   
            subtotal+=ind_result 
        j+=1
    return subtotal 

def integrand1(y,j=2,weights=weights,starts=starts,ends=ends):
         
    den=mu_multiple(y)
    length=ends[j]-starts[j]
    num=1/length
    integrand=np.log(num/den)*(weights[j]/length)      
  
    return integrand
def integrand2(y,j=1,weights=weights,starts=starts,ends=ends):
         
    den=mu_multiple(y)
    length=ends[j]-starts[j]
    num=1/length
    integrand=np.log(num/den)*(weights[j]/length)      
  
    return integrand

def integrand3(y,j=0,weights=weights,starts=starts,ends=ends):
         
    den=mu_multiple(y)
    length=ends[j]-starts[j]
    num=1/length
    integrand=np.log(num/den)*(weights[j]/length)       
  
    return integrand

def final(*arg):
    I = quad(integrand1,0.2,1.3)+quad(integrand2,0.1,1.2)+quad(integrand3,0,1) 
    I=(I[0]+I[2]+I[4])*np.log2(np.exp(1))    
    return I

print(final())