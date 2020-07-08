import os
import numpy as np
import matplotlib.pyplot as plt
if not os.path.exists("../pictures_iso"):
    	os.makedirs("../pictures_iso")
for root, dirs, files in os.walk(".", topdown=False):
	for name in dirs:
		os.chdir(name)
		
		#First we look at the hamiltonian
		Array_Plus=np.loadtxt("Av_Hamiltonian.txt")
		Array=Array_Plus[1:]  #Removes the first element
		N=16 #State the amount of residues
		H=np.zeros((N,N))   #Create an empty hamiltonian
		
		#From this created array we want to create an upper triangular matrix
		#using the indexing rule
		iu1 = np.triu_indices(N)  #Creates indices for upper triangular matrix
		H[iu1]=Array
		#and now we flip it around the diagonal axis, 
		for i in range(len(H)):
			for j in range(len(H)):
				H[j,i]=H[i,j]

		#erase the strong diagonal couplings 
		for i in range(len(H)):
			H[i,i]=0

		#Create an effective gradient
		H[0,0]=-11
		H[1,1]=15
		
		#Show the Matrix coloured 
		#This shows relative scale and intensity of present couplings
		
		fig, ax = plt.subplots()
		im = ax.imshow(H, cmap=plt.get_cmap('tab20c'), interpolation='nearest')
		fig.colorbar(im)
		ax.set_xlabel("Protein Chain Number",fontsize=16)
		ax.set_ylabel("Protein Chain Number",fontsize=16)
		plt.xticks(fontsize=12)
		plt.yticks(fontsize=12)
		plt.savefig("../../pictures_iso/hamiltonian_"+name+".png", dpi=200)
		plt.clf()
		#This first section can be applied to any of the proteins within this thesis to create the visualized Hamiltonians. The most important part is that the amount of residues and the utilized Hamiltonian text file are correct.


				
		
		#Here we plot the 1D Reponse Function
		Data=np.loadtxt('TD_Absorption.dat') #cc=0
		Data_2=np.loadtxt('TD_Absorption_2.dat') #cc=100
		plt.plot(Data[:,0],Data[:,1],c='firebrick',label='Fully Labeled Sample',linewidth=2)
		plt.plot(Data_2[:,0],Data_2[:,1],c='steelblue',label='Diluted Labeled Sample',linewidth=2)
		plt.xlabel('Time [fs]',fontsize=16)
		plt.ylabel('Response function [arb.u.]',fontsize=16)
		plt.xticks(fontsize=12)
		plt.yticks(fontsize=12)
		plt.legend()
		plt.savefig("../../pictures_iso/"+name+"_response.png", dpi=200)
		plt.clf()
		#plt.show()

		#Here we plot the 1D Absorption Spectrum
		Data=np.loadtxt('Absorption.dat')
		Data_2=np.loadtxt('Absorption_2.dat') #cc=100
		X=Data[:,0]
		Y=Data[:,1]
		X_2=Data_2[:,0]
		Y_2=Data_2[:,1]
		plt.plot(X,Y,c='firebrick',label='Fully Labeled Sample',linewidth=2)
		plt.plot(X_2,Y_2,c='steelblue',label='Diluted Labeled Sample',linewidth=2)
		# Find the highest point and the X value accociated with it
		max_Y = max(Y)  # Find the maximum y value
		max_X = X[Y.argmax()]  # Find the x value corresponding to the maximum y value
		max_Y_2=max(Y_2)
		max_X_2=X_2[Y_2.argmax()]
		print(name)
		print("Fully Labeled Sample",max_X, max_Y)
		print("Diluted Labeled Sample",max_X_2,max_Y_2)
		plt.xlabel('$\omega$ [cm$^{-1}$]',fontsize=16)
		plt.xlim(1540,1700)
		plt.ylabel('Absorption [arb.u.]',fontsize=16)
		plt.xticks(fontsize=12)
		plt.yticks(fontsize=12)
		plt.legend()
		plt.savefig("../../pictures_iso/"+name+"_absorption.png", dpi=200)
		plt.clf()
		#plt.show()
		
		
		
		
		
		os.chdir("..")
