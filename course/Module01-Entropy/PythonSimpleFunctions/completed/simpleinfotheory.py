import numpy as np
import math

"""function infocontent(p)
Computes the Shannon information content for an outcome x of a random variable
X with probability p.

Inputs:
- p - probability to compute the Shannon info content for

Outputs:
- result - Shannon info content of the probability p

Copyright (C) 2020-, Julio Correa, Joseph T. Lizier
Distributed under GNU General Public License v3
"""

def infocontent(p):
    
    # Alter the equation below to provide the correct Shannon information 
    # content:

    return -np.log2(p)

"""function entropy(p)
Computes the Shannon entropy for a probability distribution p.

Inputs:
- p - (array which much sum to 1) - a probability distribution to compute the Shannon info content for

Outputs:
- result - Shannon entropy of the probability distribution p
 
Copyright (C) 2020-, Julio Correa, Joseph T. Lizier
Distributed under GNU General Public License v3
"""
def entropy(p):  
    # First make sure the array is now a numpy array
    if type(p) != np.array:
        p = np.array(p)

    # Should we check any potential error conditions on the input?
    if (abs(np.sum(p) - 1) > 0.00001):
        raise Exception("Probability distribution must sum to 1: sum is %.4f" % np.sum(p))
    
    # We need to take the expectation value over the Shannon info content at
    # p(x) for each outcome x:
    weightedShannonInfos = p*(infocontent(p))
    # nansum ignores the nans from calling infocontent(0), but we still get the warning if an entry in p is zero
    return np.nansum(weightedShannonInfos)

#################################
# End of module 1 functions
#################################

""" function entropyempirical(xn)
Computes the Shannon entropy over all outcomes x of a random variable
X from samples x_n.

Inputs:
- xn - samples of outcomes x.
   xn is a column vector, e.g. xn = [0;0;1;0;1;0;1;1;1;0] for a binary variable.

Outputs:
- result - Shannon entropy over all outcomes
- symbols - list of unique samples
- probabilities - probabilities for each sample

Copyright (C) 2020-, Julio Correa, Joseph T. Lizier
Distributed under GNU General Public License v3
"""
def entropyempirical(xn):

    # First, error checking, and converting argument into standard form:    
    if type(xn) == list:
        xn = np.array(xn)
    if xn.ndim == 1:
        xn = np.reshape(xn,(len(xn), 1)) #reshaping our 1-dim vector to numpy format of a column vector
    [xnSamples,xnDimensions] = xn.shape
    
    # We need to work out the alphabet here.
    # The following returns a vector of the alphabet:    
    symbols = np.unique(xn, axis=0)
    
    # It would be faster to call:
    [symbols, counts] = np.unique(xn, axis=0, return_counts=True)
    
    # but we could count the samples manually below for instructive purposes:
    # Next we need to count the number of occurances of each symbol in 
    # the alphabet:
    # counts = []
    # for symbol in symbols:
    #    count = 0
    #    for row in xn:
    #        if (row==symbol).all():
    #            count += 1
    #    counts.append(count)
    # counts = np.array(counts);
    
    # Now normalise the counts into probabilities:
    probabilities = counts / xnSamples
    
    # Once we have the probabilities we can simply call our existing function:
    result = entropy(probabilities)
    
    return result, symbols, probabilities
    
""" function jointentropy(p)
Computes the joint Shannon entropy over all outcome vectors x of a vector
random variable X with probability matrix p(x) for each candidate outcome
vector x.

Inputs:
- p - probability distribution function over all outcome vectors x.
   p is a matrix over all combinations of the sub-variables of x,
where p(1,3) gives the probability of the first symbol of sub-variable
x1 co-occuring with the third symbol of sub-variable x2.
   E.g. p = [0.2, 0.3; 0.1, 0.4]. The sum over p must be 1.

Outputs:
- result - joint Shannon entropy of the probability distribution p

Copyright (C) 2020-, Julio Correa, Joseph T. Lizier
Distributed under GNU General Public License v3
"""
def jointentropy(p):
    
    # Should we check any potential error conditions on the input?

    # We need to take the expectation value over the Shannon info content at
    #  p(x) for each outcome x in the joint PDF:
    # Hint: will your code for entropy(p) work, or can you alter it slightly
    #  to make it work?
    
    joint_entropy = entropy(p)
    
    return joint_entropy

""" function jointentropyempirical(xn, yn)
Computes the Shannon entropy over all outcome vectors x of a vector random
variable X from sample vectors x_n. User can call with two such arguments 
if they don't wish to join them outside of the call.

Inputs:
- xn - matrix of samples of outcomes x. May be a 1D vector of samples
    (in which case yn is also supplied), or
    a 2D matrix, where each row is a vector sample for a multivariate X
    (in which case yn is not supplied).
- yn - as per xn, except that yn is not required to be supplied (in which
 case the entropy is only calculated over the multivariate xn variable).

Outputs:
- result - joint Shannon entropy over all samples
- symbols - list of unique joint vector samples
- probabilities - probabilities for each joint symbol

Copyright (C) 2020-, Julio Correa, Joseph T. Lizier
Distributed under GNU General Public License v3
"""
def jointentropyempirical(xn, yn=[]):
    
    # First, error checking, and converting argument into standard form:    
    xn = np.array(xn)
    # Convert to column vectors if not already:
    if xn.ndim == 1:
        xn = np.reshape(xn,(len(xn),1))
    yn = np.array(yn)
    if (yn.size > 0):
        # Convert to column vectors if not already:
        if yn.ndim == 1:
            yn = np.reshape(yn,(len(yn),1))
        [rx,cx] = xn.shape
        [ry,cy] = yn.shape
        # Check that their number of rows are the same:
        assert(rx == ry)
        # Now joint them up so we only need work with xn
        xn = np.concatenate((xn,yn), axis=1)
        
    # TRICK: Next combine the row vectors in each sample into a single 
    #  symbol (being the index from the symbols array,
    # so that we can simply compute entropy on that combined symbol
    [symbols, symbolIndexForEachSample] = np.unique(xn, axis=0, return_inverse=True)

    # And compute the entropy using our existing function:
    [result, symbols_of_indices, probabilities] = entropyempirical(symbolIndexForEachSample);

    # The order of symbols is the same as their order for the probabilities

    return result, symbols, probabilities

"""function conditionalentropy(p)

Computes the conditional Shannon entropy over all outcomes x of a random
variable X, given outcomes y of a random variable Y.
Probability matrix p(x,y) is given for each candidate outcome
(x,y).

Inputs:
- p - 2D probability distribution function over all outcomes (x,y).
   p is a matrix over all combinations of x and y,
where p(1,3) gives the probability of the first symbol of variable
x co-occuring with the third symbol of variable y.
   E.g. p = [0.2, 0.3; 0.1, 0.4]. The sum over p must be 1.

Outputs:
- result - conditional Shannon entropy of X given Y

Copyright (C) 2020-, Julio Correa, Joseph T. Lizier
Distributed under GNU General Public License v3
"""
def conditionalentropy(p):
    
    # First make sure the array is now a numpy array
    if type(p) != np.array:
        p = np.array(p)

    # Should we check any potential error conditions on the input?
    # a. Should we check p is a matrix, not a vector?
    # Actually we won't since a vector would be valid if one variable only ever took one value.
    # b. Check that the probabilities normalise to 1:
    if (abs(np.sum(p) - 1) > 0.00001):
        raise Exception("Probability distribution must sum to 1: sum is %.4f" % np.sum(p))

    # We need to compute H(X,Y) - H(X):
    # 1. joint entropy: Can we re-use existing code?
    H_XY = jointentropy(p);
    # 2. marginal entropy of Y: Can we re-use existing code?
    #  But how to get p_y???
    p_y = p.sum(axis=0); # Since y changes along the columns, summing over the x's (dimension 0 argument in the sum) will just return p(y)
    H_Y = entropy(p_y);
    
    result = H_XY - H_Y;
    return result

"""function conditionalentropyempirical(xn, yn)
Computes the conditional Shannon entropy over all samples xn of a random
variable X, given samples yn of a random variable Y.

Inputs:
- xn - matrix of samples of outcomes x. May be a 1D vector of samples, or
    a 2D matrix, where each row is a vector sample for a multivariate X.
- yn - matrix of samples of outcomes x. May be a 1D vector of samples, or
    a 2D matrix, where each row is a vector sample for a multivariate Y.
    Must have the same number of rows as X.

Outputs:
- result - conditional Shannon entropy of X given Y

Copyright (C) 2020-, Julio Correa, Joseph T. Lizier
Distributed under GNU General Public License v3
"""
def conditionalentropyempirical(xn, yn):
    
    # First, error checking, and converting argument into standard form:    
    xn = np.array(xn)
    # Convert to column vectors if not already:
    if xn.ndim == 1:
        xn = np.reshape(xn,(len(xn),1))
    yn = np.array(yn)
    if yn.ndim == 1:
        yn = np.reshape(yn,(len(yn),1))
    [rx,cx] = xn.shape
    [ry,cy] = yn.shape

    # Should we check any potential error conditions on the input?
    # Check that their number of rows are the same:
    assert(rx == ry)
    
    # We need to compute H(X,Y) - H(X):
    # 1. joint entropy: Can we re-use existing code?
    (H_XY, xySymbols, xyProbs) = jointentropyempirical(xn, yn);
    # 2. marginal entropy of Y: Can we re-use existing code?
    (H_Y, ySymbols, yProbs) = entropyempirical(yn);
    
    result = H_XY - H_Y;
    return result

#################################
# End of module 2 functions
#################################

"""function mutualinformation(p)
Computes the mutual information over all outcomes x of a random
variable X with outcomes y of a random variable Y.
Probability matrix p(x,y) is given for each candidate outcome
(x,y).

Inputs:
- p - 2D probability distribution function over all outcomes (x,y).
   p is a matrix over all combinations of x and y,
where p(1,3) gives the probability of the first symbol of variable
x co-occuring with the third symbol of variable y.
   E.g. p = [0.2, 0.3; 0.1, 0.4]. The sum over p must be 1.

Outputs:
- result - mutual information of X with Y

Copyright (C) 2020-, Julio Correa, Joseph T. Lizier
Distributed under GNU General Public License v3
"""
def mutualinformation(p):
    
    # First make sure the array is now a numpy array
    if type(p) != np.array:
        p = np.array(p)

    # Should we check any potential error conditions on the input?
    # a. Should we check p is a matrix, not a vector?
    # Actually we won't since a vector would be valid if one variable only ever took one value.
    # b. Check that the probabilities normalise to 1:
    if (abs(np.sum(p) - 1) > 0.00001):
        raise Exception("Probability distribution must sum to 1: sum is %.4f" % np.sum(p))

    # We need to compute H(X) + H(Y) - H(X,Y):
    # 1. joint entropy:
    H_XY = jointentropy(p)

    # 2. marginal entropy of X:
    # But how to get p_x???
    p_x = p.sum(axis=1); # Since x changes along the rows, summing over the y's (dimension 1 argument in the sum) will just return p(x)
    H_X = entropy(p_x);

    # 2. marginal entropy of Y:
    # But how to get p_y???
    p_y = p.sum(axis=0); # Since y changes along the columns, summing over the x's (dimension 0 argument in the sum) will just return p(y)
    H_Y = entropy(p_y);

    result = H_X + H_Y - H_XY
    
    return result

"""function mutualinformationempirical(xn,yn)
Computes the mutual information over all samples xn of a random
variable X with samples yn of a random variable Y.

Inputs:
- xn - matrix of samples of outcomes x. May be a 1D vector of samples, or
    a 2D matrix, where each row is a vector sample for a multivariate X.
- yn - matrix of samples of outcomes x. May be a 1D vector of samples, or
    a 2D matrix, where each row is a vector sample for a multivariate Y.
    Must have the same number of rows as X.

Outputs:
- result - mutual information of X with Y
- xySymbols - list of unique joint vector samples
- xyProbs - probabilities for each joint symbol
- xSymbols - list of unique x samples
- xProbs - probabilities for each x symbol
- ySymbols - list of unique y samples
- yProbs - probabilities for y symbol

Copyright (C) 2020-, Julio Correa, Joseph T. Lizier
Distributed under GNU General Public License v3
"""
def mutualinformationempirical(xn,yn):
    
    # First, error checking, and converting argument into standard form:    
    xn = np.array(xn)
    # Convert to column vectors if not already:
    if xn.ndim == 1:
        xn = np.reshape(xn,(len(xn),1))
    yn = np.array(yn)
    if yn.ndim == 1:
        yn = np.reshape(yn,(len(yn),1))
    [rx,cx] = xn.shape
    [ry,cy] = yn.shape

    # Should we check any potential error conditions on the input?
    # Check that their number of rows are the same:
    assert(rx == ry)

    # We need to compute H(X) + H(Y) - H(X,Y):
    # 1. joint entropy:
    (H_XY, xySymbols, xyProbs) = jointentropyempirical(xn, yn); # How to compute this empirically ...?
    # 2. marginal entropy of Y: (call 'joint' in case yn is multivariate)
    (H_Y, ySymbols, yProbs) = jointentropyempirical(yn)
    # 3. marginal entropy of X: (call 'joint' in case yn is multivariate)
    (H_X, xSymbols, xProbs) = jointentropyempirical(xn);
    
    result = H_X + H_Y - H_XY;
    return result, xySymbols, xyProbs, xSymbols, xProbs, ySymbols, yProbs

#################################
# End of module 3 functions
#################################

"""function conditionalmutualinformation(p)
Computes the mutual information over all outcomes x of a random
variable X with outcomes y of a random variable Y, conditioning on 
outcomes z of a random variable Z.
Probability matrix p(x,y,z) is given for each candidate outcome
(x,y,z).

Inputs:
- p - 3D probability distribution function over all outcomes (x,y,z).
   p is a matrix over all combinations of x and y and z,
where p(0,2,1) gives the probability of the first symbol of variable
x co-occuring with the third symbol of variable y and the second
symbol of z.
   The sum over p must be 1.
   E.g.:
     p[0,:,:] = [[0.114286, 0.171429], [0.057143, 0.228571]]
     p[1,:,:] = [[0.171429, 0.114286], [0.028571, 0.114286]]

Outputs:
- result - conditional mutual information of X with Y given Z
 
Copyright (C) 2020-, Julio Correa, Joseph T. Lizier
Distributed under GNU General Public License v3
"""
def conditionalmutualinformation(p):
    
    # First make sure the array is now a numpy array
    if type(p) != np.array:
        p = np.array(p)

    # Should we check any potential error conditions on the input?
    # a. Check that we have 3 dimensions. We allowed one dimension to be null
    #  for MI and conditional entropy, but for CMI we can't tell which is missing
    if (p.ndim != 3):
        raise Exception("Probability distribution must have 3 dimensions for CMI")
    # b. Check that the probabilities normalise to 1:
    if (abs(np.sum(p) - 1) > 0.00001):
        raise Exception("Probability distribution must sum to 1: sum is %.4f" % np.sum(p))

    # We need to compute H(X|Z) + H(Y|Z) - H(X,Y|Z).
    # But our conditional entropy calculator won't do H(X,Y|Z) since it doesn't accept a joint probability for X,Y.
    # So, easier to rewrite as:
    #  H(X,Z) - H(Z) + H(Y,Z) - H(Z) - H(X,Y,Z) + H(Z)
    #  = H(X,Z) - H(Z) + H(Y,Z) - H(X,Y,Z)
    
    # 1. joint entropy:
    H_XYZ = jointentropy(p)
    
    # 2. entropy of X,Z:
    #  But how to get p_xz???
    # Sum p over the y's (2nd dimension argument in the sum) will just return p(x,z) terms.
    p_xz = p.sum(axis=1)
    H_XZ = jointentropy(p_xz)

    # 3. entropy of Y,Z:
    #  But how to get p_yz???
    # Sum p over the x's (1st dimension argument in the sum) will just return p(y,z) terms.
    p_yz = p.sum(axis=0)
    H_YZ = jointentropy(p_yz)

    # 4. marginal entropy of Z:
    #  But how to get p_z???
    # Sum p_xz over the x's (1st dimension argument in the sum) will just return p(z) terms.
    p_z = p_xz.sum(axis=0)
    H_Z = jointentropy(p_z)
    
    result = H_XZ - H_Z + H_YZ - H_XYZ
    return result

"""function conditionalmutualinformationempirical(xn,yn,zn)
Computes the mutual information over all samples xn of a random
variable X with samples yn of a random variable Y, conditioning on 
samples zn of a random variable Z.

Inputs:
- xn - matrix of samples of outcomes x. May be a 1D vector of samples, or
    a 2D matrix, where each row is a vector sample for a multivariate X.
- yn - matrix of samples of outcomes y. May be a 1D vector of samples, or
    a 2D matrix, where each row is a vector sample for a multivariate Y.
    Must have the same number of rows as X.
- zn - matrix of samples of outcomes z. May be a 1D vector of samples, or
    a 2D matrix, where each row is a vector sample for a multivariate Z
    which will be conditioned on.
    Must have the same number of rows as X.

Outputs:
- result - conditional mutual information of X with Y, given Z

Copyright (C) 2020-, Julio Correa, Joseph T. Lizier
Distributed under GNU General Public License v3
"""
def conditionalmutualinformationempirical(xn, yn, zn):
    
    # First, error checking, and converting argument into standard form:    
    xn = np.array(xn)
    # Convert to column vectors if not already:
    if xn.ndim == 1:
        xn = np.reshape(xn,(len(xn),1))
    yn = np.array(yn)
    if yn.ndim == 1:
        yn = np.reshape(yn,(len(yn),1))
    zn = np.array(zn)
    if zn.ndim == 1:
        zn = np.reshape(zn,(len(zn),1))
    [rx,cx] = xn.shape
    [ry,cy] = yn.shape
    [rz,cz] = zn.shape

    # Should we check any potential error conditions on the input?
    # Check that their number of rows are the same:
    assert(rx == ry)
    assert(rx == rz)

    # We need to compute H(X|Z) + H(Y|Z) - H(X,Y|Z):
    # 1. conditional joint entropy:
    H_XY_given_Z = conditionalentropyempirical(np.append(xn, yn, axis=1),zn); # How to compute this empirically ...?
    # 2. conditional entropy of Y:
    H_Y_given_Z = conditionalentropyempirical(yn,zn) # How to compute this empirically ...?
    # 3. conditional entropy of X:
    H_X_given_Z = conditionalentropyempirical(xn,zn) # How to compute this empirically ...?
    
    # Alternatively, note that we could compute I(X;Y,Z) - I(X;Z)
    
    result = H_X_given_Z + H_Y_given_Z - H_XY_given_Z;
    return result

