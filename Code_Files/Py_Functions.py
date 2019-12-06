from numpy import array, zeros, log10
import matplotlib.pyplot as plt
from timeit import default_timer as timer

def readarrays(filename):
	#start = timer()
	values = open(filename, "r")
	#print values.read()	
	lines = values.readlines()
	#end = timer()
	#print (end-start)


	#Counting
	C = 0
	D = 0
	Dims = [] 
	A = []

	#start = timer()
	for i in lines:
		if i != "\n":
			D += 1
		if i == "\n":
			C += 1
			Dims.append(D)
			A.append(zeros(D))
			D = 0
	#end = timer()
	#print (end-start)

	#start = timer()	
	#Filling
	F = 0
	G = 0
	for i in lines:
		if i != "\n":
			A[F][G] = i
			G += 1
		if i == "\n":
			F += 1
			G = 0 
	#end = timer()
	#print (end-start)
	values.close()
	return A,len(A)

def readmatrices(filename):
	#start = timer()
	values = open(filename, "r")
	#print values.read()	
	lines = values.readlines()
	#end = timer()
	#print (end-start)


	#Counting
	C = 0
	D = 0
	Dims = [] 
	A = []

	#start = timer()
	
	for i in range(len(lines)):
		if lines[i] != "\n":
			D += 1
		if lines[i] == "\n":
			C += 1
			Dims.append(D)
			A.append(zeros(shape=(D,len(lines[i-1].split()))))
			D = 0	

	
	#start = timer()	
	#Filling
	F = 0
	G = 0
	for i in lines:
		if i != "\n":
			for j in range(len(i.split())):			
				A[F][G][j] = i.split()[j]
			G += 1
		if i == "\n":
			F += 1
			G = 0 
	#end = timer()
	#print (end-start)
	values.close()
		
	return A,len(A)

def plot_median_d(filename,save=False):
	median = readarrays(filename)[0][0]
	D = filename.split("_") ; D = int(D[2])
	filename = filename.split("/") ; filename = filename[2].split(".txt") ; filename = filename[0]
	N = log10( (10**(D-2)) * array(range(1,len(median)+1)) )

	plt.figure()
	plt.title(filename)
	plt.plot(N,median,'.',N,median,'r')
	plt.xlabel("$\log_{10}(MC)$")
	plt.ylabel("$\mu_{1/2}$")
	if save:
		Folder = "Plots_d_D_%s/" % D
		plt.savefig("../Plots/Plots_d/" + Folder + filename + ".png")	

	plt.close()


