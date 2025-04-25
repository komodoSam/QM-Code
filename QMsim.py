# Quantum Mechanics Simulation
import sys
import numpy as np
import math
import scipy.sparse as sp
import matplotlib.pyplot as plt
from numpy import linalg as la

# constants
pi = np.pi

Nx = 1000  # made this five just so I can see the full matrix (original is 1,000)
#amount of segments that graph split into 

dx = 2*pi/Nx  
k = 40  

def V(x):
    return 0.5*k*x**2   # Potential energy for k = 40 (I defined the variable)

# matrices
mat = np.zeros((Nx+1, Nx+1))   # maybe a square matrix (of zeroes probably) with Nx+1 rows and columns?

for i in range(Nx+1):    # nested for loop: 2D array / matrix traversal
    for j in range(Nx+1): 
        if i == j: mat[i,j] = -2
        if i == j + 1: mat[i,j] = 1   # the top left to bottom right diagonal made of "-2"s
        if i == j - 1: mat[i,j] = 1   # and the above and below diagonals are comprised of "1"s  (rest are zeroes)

lap = (-1/dx**2)*mat/20  # makes another matrix of same size and with the same diagonals, but different values
#print(lap)
# has negative zeroes?? maybe just syntax thing
# Stands for Laplacian: Laplacian - derivative^2

pot = [V((2*pi*i/Nx)-pi) for i in range(Nx+1)]  # calculating the energy of each position on matrix
#print(pot)

pot_mat = np.zeros((Nx+1, Nx+1))  # potential matrix like mat (same size)

for i in range(Nx+1):
    for j in range(Nx+1):
        if i == j: pot_mat[i,j] = pot[i]   # getting flashbacks, but this one uses pot's values for the diagonal

ham = lap + pot_mat  # hamiltonian - operator (matrix): K and U correspond to lap and pot
#Physically its just the matrix sum of lap and pot_mat
#print(ham)

w,v = la.eig(ham)  # eigenvalues or eigenvectors idk 
#print(w)
#print(v)

idx = w.argsort()[::1]
#stores eigen values and eigen vectors in relation to another (lowest to highest energy)
eigenValues = w[idx]    # knew it :D
eigenVectors = v[:,idx] 

#idx -> array([ 21, 167, 131, ..., 1, 0, 181], dtype=int64)  # labelled as Out[13] on paper -> output probably



x1 = -390  #if they were same as A max, then the denominator would go to 0    
x2 = 390  
# From what I can tell these are the range for the points graphed on the plot
#I altered them to be able to get rid of "-500"
Amax = 400  #Max Amplitude 

#class_prob = [0.65/np.sqrt(400**2-(i-500)**2) for i  in range(x1, x2)]   #old classical prob values
class_prob = [(1/2*pi) * (1/np.sqrt(Amax**2 - i**2)) * 100 for i in range(x1, x2)]

#plt.plot(range(x1, x2), class_prob)    # the first graph displayed in the first plot
#plt.plot(eigenVectors[:,1])    # represents the second graph shown in the first plot
#plt.show()
# the second graph and the two below this are the same but with different numbers of oscillations

#plt.plot(eigenVectors[:,50])  #This isn't actually part of the code, I just modified the other graph
#plt.show()
#This makes a plot similar to the second picture 

# plotting different eigenvectors 
# colon represents all possible values for the number

#plt.plot(eigenVectors[:,150]**2)  # Maybe models wavefunction -> Quantum walk graph, would explain why it's squared
#plt.show()
#Note: In QM, the eigenvalues of the Hamiltonian matrix correspond to the energy levels that are observable in experiments
#the eigenvectors correspond to the wave functions. 


# Output: [<matplotlib.lines.Line2D at 0x18d0029da00>]

xvals = np.zeros(Nx+1) # plotting the energies
spec = np.sort(w)  # plotting some of the lowest energies (energies evenly spaced conveniently)
#plt.ylim(0,10)
#plt.scatter(xvals, spec)  # don't know what this scatter plot represents D:
#plt.show()

diff = [spec[n+1] - spec[n] for n in range(Nx)]
diff2 = [diff[n+1] - diff[n] for n in range(Nx - 1)]   #don't know what these things represent either

xvals = np.zeros(len(diff2))
spec = np.sort
#plt.plot(diff2[:200])
#plt.show()

#print(v[[0]])


def V_2(x):
    return 50*x**2 + 50*(1+np.cos(2*pi*x))   #Why do we need a second potential function?
# method for any potential, allows for more complicated potentials

Nx = 1000
dx = 2*pi/Nx   # redefining old variables

phi_rescaled = [2*pi*i/Nx-pi for i in range(Nx + 1)]
phi_vals = [2*pi*i/Nx-pi for i in range(Nx + 1)]   # represents the position the particle could be  (Phase of the particle)
V2_vals = [V_2(phi) for phi in phi_vals]
plt.plot(phi_vals, V2_vals)
plt.ylim(-50, 200)
plt.show()

mat = np.zeros((Nx+1, Nx+1))

for i in range(Nx+1):
    for j in range(Nx+1):
        if i == j: pot_mat[i,j] = pot[i]

ham = lap + pot_mat
w,v = la.eig(ham)

idx = w.argsort()[::1]
eigenValues = w[idx]
eigenVectors = v[:,idx]

V2_rescaled = [V2 for V2 in V2_vals]

#plt.plot(range(Nx+1), V2_rescaled)
#for i in range(0, 2):
    #plt.plot(75*eigenVectors[:,i]+eigenValues[i])
#plt.show()