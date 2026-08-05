from scipy.integrate import quad
import numpy as np

mu=[0.4,0.5,0.8] #Default values are from Ross paper
si=[0.2, 0.3, 0.25]
weights=[0.3,0.5,0.2]
def exp(y,m,s):
    frac=((y-m)**2)/(s**2)
    exp=np.exp(-0.5*frac)/s
    exp=exp/np.sqrt(2*np.pi)
    return exp#

def denom(y,weights=weights,list_of_mu=mu,list_of_sigma=si):
    i=0
    e_values=np.zeros(len(weights))
    for i in range(len(weights)):
        e_values[i]=exp(y,list_of_mu[i],list_of_sigma[i])*weights[i] 
    return np.sum(e_values)       
       
#print(type(denom(1)))
def integrand(y, weights=weights,list_of_mu=mu,list_of_sigma=si):
    j=0
    integrand=0
    while j<len(weights):
        num=exp(y,list_of_mu[j], list_of_sigma[j])*weights[j]
        den=denom(y)*weights[j]
        contribution=np.log(num/den)*num
        integrand+= contribution 
        j=j+1    
    return integrand
#def integrand_test(y,m=mu[2],s=si[2]):
#    return exp(y,m,s)      
#J=quad(integrand_test,-5,5)
I = quad(integrand, -5,5)
I=I[0]*np.log2(np.exp(1)) 
print(str(I)+' bits')
