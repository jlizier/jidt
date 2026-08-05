#0.Setup____________________________
#0.0. import libraries and setup paths
import sys
from jpype import *
import numpy as np 
import matplotlib.pyplot as plt
import DataGenerators as dg
import integral_gaussian as ig
import integral_square as igs 

sys.path.append("C:/Users/gmdf9/Dropbox/CSYS AND MATHS/Capstone JL/demos/python")
addClassPath("C:/Users/gmdf9/Dropbox/CSYS AND MATHS/Capstone JL/infodynamics.jar")

print('What type of output would you like to use?')
print('0: sends results to a .csv file in your working directory')
print('1: produces pyplot with error bars')
print('2: produces pyplot with shaded error area')
type_output=int(input(''))
print('What kind of data would you like to model?')
print('0: Gaussian data')
print('1: Square wave data')
print('2: Correlated Gaussian Data (CMI ONLY)')
datatype=int(input(''))
if datatype==2:
    calctype=1
else:
    print('Which calculation are you performing?')
    print('0: Mutual information')
    print('1: Conditional Mutual Information')
    calctype=int(input(''))



if not isJVMStarted(): #important; otherwise need to restart IDE to run new computation
    startJVM("-ea", convertStrings=False)  

# 0.1. Construct the calculator:
if calctype==0:
    calcClass = JPackage("infodynamics.measures.mixed.kraskov").MutualInfoCalculatorMultiVariateWithDiscreteKraskov
    calc = calcClass()
else:
    calcClass = JPackage("infodynamics.measures.mixed.kraskov").ConditionalMutualInfoCalculatorMultiVariateWithDiscreteKraskov
    calc = calcClass()
    
    

#0.2. Set up lists and the maximum k=value we are interested in, for eventual plot.
L=10
if type_output==0:
    final_results=[]
elif type_output==1:    
    x_axis=np.zeros(L)
    y_axis=np.zeros(L)
    y_errs=np.zeros(L)
else:
    x_axis=np.zeros(L)
    y_axis=np.zeros(L)
    error_array=np.zeros((2,L))

#0.3. Setup data collection structures before generating data:
reps=100
size=400 
list_abund=[2,10,5] # determines how many discrete values & their relative prevalence
ab_norm=list_abund/sum(np.array(list_abund)) #normalizing
#0.3.0 Set up for square wave data

#0.3.1 Set up for Gaussian data. Default values are from Ross paper




#2. MAIN METHOD IN JIDT(stages 2.1-2.4 are repeated)___________________________________________________________-
l=0
while l<int(L):
    # 2.0 Set any properties to non-default values:
    l=l+1 #Normally would increment count at the end, but we REALLY do not want k=0
    calc.setProperty("k",str(l));
    results=[]
    j=0
    while j<reps:
        if datatype==0:
            (x,y,r)=dg.GenerateSingleGaussianData(ab_norm=ab_norm,size=size) # by default, uses the lists from 1b.
            
        elif datatype==1:
            (x,y,r)=dg.GenerateSquareWaveData(ab_norm=ab_norm,size=size)
            
        else:
            (x,y,z,r)=dg.GenerateCorrelatedGaussianData(ab_norm=ab_norm)
            
        if calctype==0:
            
            dsc= JArray(JInt, 1)(x) #VERY VERY important not to use numpy arrays
            cts = JArray(JDouble, 1)(y) #Also important to distinguish JInt from JDouble      
            calc.initialise(1,r)            
            calc.setObservations(cts,dsc) 
        else:
            cond=JArray(JDouble,1)(z)
            dsc= JArray(JInt, 1)(x)
            cts = JArray(JDouble, 1)(y)
            calc.initialise(1,r,1)       
       
            calc.setObservations(cts,dsc,cond) 
        
        # 2.3. Compute the estimate for each rep:
        result = calc.computeAverageLocalOfObservations()
        results.append(result)
        j=j+1
        
    #2.4.Amalgamate the reps into usable form for each k
    results=np.array(results)*np.log2(np.exp(1))    
    bits=np.average(results)
    sd=np.std(results)
    print(bits,sd)

    if type_output==0:
            final_results.append((l,bits, sd))
    elif type_output==1:   
            x_axis[l-1]=l
            y_axis[l-1]=bits
            y_errs[l-1]=2*sd
            
    else:
        x_axis[l-1]=l
        y_axis[l-1]=bits
        (perc10,perc90)=np.percentile(results,(10,90))
        error_array[0][l-1]=perc10
        error_array[1][l-1]=perc90
#3. Output______________________________________________
if datatype==0:
    ig_to_plot=ig.final()
elif datatype==1:
    ig_to_plot=igs.final()

if type_output==0:
    np.savetxt('output.csv', final_results, delimiter='\t')
elif type_output==1:
    
    plt.xlabel('k')
    plt.ylabel('MI (bits)')
    plt.axis([0.9,L,0,0.4])
    plt.axhline(y=ig_to_plot, color='blue', linestyle='--')
    plt.errorbar(x_axis, y_axis,yerr=y_errs,ecolor='black')
    plt.plot()
    
elif type_output==2:
    plt.xlabel('k')
    plt.ylabel('MI (bits)')
    plt.axis([0.9,L,0,0.4])
    plt.axhline(y=ig_to_plot, color='blue', linestyle='--')
    plt.fill_between(x_axis, error_array[0], error_array[1],color='gray')
    plt.plot(x_axis,y_axis,color='black')

    



