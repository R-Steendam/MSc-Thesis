import numpy as np
import matplotlib.pyplot as plt
from vpython import *

# This subroutine visualizes a 3D structure using vpython
# The vpython library must be available

# Visualize the bare structure with white arrows for the transition dipoles
def visual(x,mu,N,scale):
    for i in range(N):
      arrow(pos=vector(x[i,0],x[i,1],x[i,2]),axis=vector(scale*mu[i,0],scale*mu[i,1],scale*mu[i,2]),)

# Visualize a selected exciton state with white arrows showing the phase and magnitude of the wave function. The molecules are shown as small spheres.
def visual_exciton(x,mu,c,index,N,scale):
    for i in range(N):
        # Show molecule
        sphere(pos=vector(x[i,0],x[i,1],x[i,2]),radius=0.5)
        # Show associated transition dipole
        arrow(pos=vector(x[i,0],x[i,1],x[i,2]),axis=vector(c[i,index]*scale*mu[i,0],c[i,index]*scale*mu[i,1],c[i,index]*scale*mu[i,2]))
        
       
# This subroutine calculates stick spectra and spectra convoluted with
# a Gaussian. Provide the Hamiltonian, H, the transition-dipoles, mu,
# the size of the Hamiltonian, N, and the standard deviation for the
# Gaussian convolution, sigma. Plots are generated. This subroutine is
# intended as a way to test the Hamiltonians generated in this package.


dip = np.loadtxt("S2PNE_Dipole.txt") #Load the copmlete Dipole Files
# dip = mu_plus[1:] #remove the first integer
N = 84 #The amount of C=O couplings within the protein, including sidechains
sigma= 1 #can be varied.

#First let's create the necessary matrix H
Array_Plus=np.loadtxt("Av_S2pne.txt") #Load the AmideIMaps Energy files
Array=Array_Plus[1:]  #Removes the first element

H=np.zeros((N,N))     #Create emtpy Hamiltonian

#From this created array we want to create an upper triangular matrix
#using the indexing rule
iu1 = np.triu_indices(N)  #Creates indices for upper triangular matrix
H[iu1]=Array
#and now we flip it around the diagonal axis, 
for i in range(len(H)):
    #H[i,i]=1650
    for j in range(len(H)):
        H[j,i]=H[i,j]


def absorption(H,N,sigma):
  # Diagonalize Hamiltonian
  # E is the eigenvalue of the Hamiltonian
  # c is NxN matrix containing the eigenvectors
  E,c=np.linalg.eigh(H)
  # Make spectrum
  bins=1000
  EminA=np.min(E)
  EmaxA=np.max(E)
  # Ensure that the resolution is bigger than the witdh of each peak
  if (EmaxA-EminA)/sigma>1000:
      bins=int(np.floor((EmaxA-EminA)/sigma))
  # Add 10% of full bandwidth and three times the convolution width on each side
  Emin=EminA-0.1*(EmaxA-EminA)-3*sigma
  Emax=EmaxA+0.1*(EmaxA-EminA)+3*sigma
  dE=(Emax-Emin)/bins
  Ex=np.linspace(Emin,Emax,bins)
  #Create emtpy bins for the data to add into
  Ey=np.zeros(bins)
  mu=np.zeros((N,3))
  mu[:,0]=dip[1:N+1]
  mu[:,1]=dip[N+1:2*N+1]
  mu[:,2]=dip[2*N+1:3*N+1]
  for n in range(N):
    bin=int(round((E[n]-Emin)/dE))
    Emu=np.zeros(3)
    Emu[0]=np.inner(c[:,n],dip[1:N+1])
    Emu[1]=np.inner(c[:,n],dip[N+1:2*N+1])
    Emu[2]=np.inner(c[:,n],dip[2*N+1:3*N+1])

    #The y value is normalized Any vector, when normalized, 
    #only changes its magnitude, not its direction
    Ey[bin]=Ey[bin]+np.linalg.norm(Emu)**2
    if np.linalg.norm(Emu)**2>20:
        print('Eigen State '+str(n)+' has Emu '+str(np.linalg.norm(Emu)**2))
        print('Eigen Frequency '+str(E[n])+' cm-1')
        print(Emu[0],"&",Emu[1],"&",Emu[2])
    if np.linalg.norm(Emu)**2>10 and E[n]>1660:
        print('Eigen State '+str(n)+' has Emu '+str(np.linalg.norm(Emu)**2))
        print('Eigen Frequency '+str(E[n])+' cm-1')
        print(Emu[0],"&",Emu[1],"&",Emu[2])
    #print(Emu)
  
  plt.plot(Ex,Ey,c='royalblue')
  plt.xlabel("Frequency",fontsize=16)
  plt.ylabel("Absorption [arb.u.]",fontsize=16)
  plt.xticks(fontsize=16)
  plt.yticks(fontsize=16)
  #plt.show()

  # Convolute spectrum
  # This is done to remove Noise from the data and create a smoother result
  # First create normalized Gaussian centered in the middle
  Cx=Ex-(Emax+Emin)/2 # Make new axis with value zero in the middle of array
  Cy=np.exp(-Cx**2/2/sigma**2)/np.sqrt(2*np.pi*sigma**2)
  plt.plot(Cx,Cy,c='dimgrey')
  plt.xlabel("Frequency",fontsize=16)
  plt.ylabel("Absorption [arb.u.]",fontsize=16)
  plt.xticks(fontsize=14)
  plt.yticks(fontsize=14)
  #plt.show()

  # Do the actual convolusion
  Ny=np.convolve(Ey,Cy,mode='same')
  plt.plot(Ex,Ny,c='firebrick')
  plt.xlabel("Frequency",fontsize=16)
  plt.ylabel("Absorption [arb.u.]",fontsize=16)
  plt.xticks(fontsize=14)
  plt.yticks(fontsize=14)
  #plt.show()

  # Plot everything in one final plot
  plt.plot(Ex,Ey/np.max(Ey),label='Absorption spectra',c='royalblue')
  plt.plot(Ex,Cy/np.max(Cy),label='Gaussian used for convolution',c='dimgrey')
  plt.plot(Ex,Ny/np.max(Ny),label='Convoluted spectra',c='firebrick')
  plt.xlabel("Frequency",fontsize=16)
  plt.ylabel("Absorption [arb.u.]",fontsize=16)
  plt.xticks(fontsize=16)
  plt.yticks(fontsize=16)
  plt.legend()
 # plt.show()
  return c,mu

x=np.loadtxt("position.txt")

c,mu=absorption(H,N,sigma)
#visual_exciton(x,mu,c,79,N,10) #Belangrijkste toestanden
##Probeer de twee pijlen te combineren met een ander kleurtje
visual(x,mu,9,3)        
