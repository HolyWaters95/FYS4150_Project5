from Py_Functions import readarrays, plot_median_d, plot_prob_distribution_d, log_arrays, nlog_arrays
import Py_Functions as pf
from numpy import array, zeros, linspace, log, log10, exp, sort, polyfit, polyval
import matplotlib.pyplot as plt
import matplotlib.image as pmg

import scipy.stats as stats
'''
N = range(1,501)
Money = []

#filename = "Money_distributions_savings_L_0.250000"
filename = "Test_sampling"

for R in range(4):
	M_R = readarrays("../build-Laptop_Project-Desktop_Qt_5_13_0_MinGW_64_bit-Debug/%s_%s.txt" % (filename,R))[0]
	for i in range(0,len(M_R)):
		Money += list(M_R[i])


w, m_intervals, patches = plt.hist(Money,len(N),density=True)

'''
'''
coeff = polyfit(m_intervals[:-1],w,10)
W = polyval(coeff,m_intervals[:-1])
plt.plot(m_intervals[:-1],W,'r')

fit_alpha, fit_loc, fit_beta=stats.gamma.fit(w)

plt.plot()
'''

'''
#plt.plot(m_intervals[:-1],w,'r')

#plt.savefig("histogram_a.png")

lw = zeros(len(w))
for i in range(len(lw)):
	if w[i] != 0:	
		lw[i] = log(w[i])
	else:
		lw[i] = w[i]

plt.figure()
plt.plot(log10(m_intervals[:-1]),log10(w))
#plt.savefig("log_w_a.png")

slope = (lw[100]-lw[1])/(m_intervals[100]-lw[1])
print slope

plt.show()


'''

D = 3


filenames = []
filestart = "../Results/Median_D_%s" % D

Nvalues = ["500", "1000"]
Lvalues = ["0.250000", "0.500000", "0.900000"]
avalues = ["0.500000", "1.000000", "1.500000", "2.000000"]
gvalues = ["1.000000", "2.000000", "3.000000", "4.000000"]
'''
#task a
filenames.append(filestart + "_N_" + Nvalues[0] + ".txt")

#task c
for L in Lvalues:
	filenames.append(filestart + "_N_" + Nvalues[0] + "_L_" + L + ".txt")
'''

'''
#task e
for L in ["0.000000","0.500000"]:
	for a in ["1.000000","2.000000"]:
		for g in gvalues:
			filenames.append(filestart + "_N_" + Nvalues[1] + "_L_" + L + "_a_" + a + "_g_" + g + ".txt")

'''



yn = raw_input("Do you want to plot medians for task d)? y/n \n")
if yn == "y":
	filenames = []
	D = raw_input("Exponent of D? \n")
	filestart = "../Results/Median_D_%s" % D	
	for N in Nvalues:
		for L in ["0.000000","0.500000"]:
			for a in avalues:
				filenames.append(filestart + "_N_" + N + "_L_" + L + "_a_" + a + ".txt")

	for f in filenames:
		plot_median_d(f,save=False)

yn = raw_input("Do you want to plot probability distributions for d)? y/n \n")
if yn == "y":
	filenames = []
	D = raw_input("Exponent of D? \n")
	filestart = "../Results/Test_Money_distributions_D_%s" % D	
	for N in Nvalues:
		for L in ["0.000000","0.500000"]:
			for a in avalues:
				filenames.append(filestart + "_N_" + N + "_L_" + L + "_a_" + a + ".txt")
				#if N == "1000" and L == "0.000000" and a == "2.000000":
				#	filenames.append(filestart + "_N_" + N + "_L_" + L + "_a_" + "9.000000" + ".txt")	
	

	W = [] ; M = []
	for i in range(len(filenames)):
		f = filenames[i]
		w, m = plot_prob_distribution_d(f,save=False)
		W.append(w) ; M.append(m[:-1])
		
	i = 0	
	while i < len(W):
		
		pf.plot_loglogW_d(M[i:i+4],W[i:i+4],filenames[i:i+4],save=False)
		pf.Pareto_dist_d(M[i:i+4],W[i:i+4],filenames[i:i+4],save=False)
		pf.Gibbs_dist_d(M[i:i+4],W[i:i+4],filenames[i:i+4],save=False)		
		i += 4

yn = raw_input("Do you want to plot medians for task e)? y/n \n")
if yn == "y":
	filenames = []
	D = raw_input("Exponent of D? \n")
	filestart = "../Results/Median_D_%s_N_1000" % D	
	for L in ["0.000000","0.500000"]:
		for a in ["1.000000","2.000000"]:
			for g in gvalues:
				filenames.append(filestart + "_L_" + L + "_a_" + a + "_g_" + g + ".txt")
	filenames.append("../Results/Median_D_8_N_1000_L_0.500000_a_2.000000_g_1.000000.txt")
	#for f in filenames:
	#	pf.plot_median_e(f,save=True)
	pf.plot_median_e(filenames[-1],save=True)

'''
filename = "../Results/Money_distributions_D_7_N_1000_L_0.000000_a_2.000000.txt"
D, N, L, a = pf.extract_parametres_d(filename)
print D, N, L, a
'''

yn = raw_input("Do you want to plot probability distributions for e)? y/n \n")
if yn == "y":
	filenames = []
	D = raw_input("Exponent of D? \n")
	filestart = "../Results/TEST_Money_distributions_D_%s_N_1000" % D	
	for L in ["0.000000","0.500000"]:
		for a in ["1.000000","2.000000"]:
			for g in gvalues:
				filenames.append(filestart + "_L_" + L + "_a_" + a + "_g_" + g + ".txt")
				#if N == "1000" and L == "0.000000" and a == "2.000000":
				#	filenames.append(filestart + "_N_" + N + "_L_" + L + "_a_" + "9.000000" + ".txt")	
	

	W = [] ; M = []
	for i in range(len(filenames)):
		f = filenames[i]
		w, m = pf.plot_prob_distribution_e(f,save=False)
		W.append(w) ; M.append(m[:-1])
		
	i = 0	
	while i < len(W):
		
		pf.plot_loglogW_e(M[i:i+4],W[i:i+4],filenames[i:i+4],save=False)
		pf.Pareto_dist_e(M[i:i+4],W[i:i+4],filenames[i:i+4],save=True)
		#pf.Gibbs_dist_e(M[i:i+4],W[i:i+4],filenames[i:i+4],save=False)		
		i += 4




'''
w,m = plot_prob_distribution_d("../Results/Money_distributions_D_7_N_1000_L_0.000000_a_2.000000_g_2.000000.txt")

logm, logw = log_arrays(m[:-1],w)
plt.plot(logm,logw)
plt.show()
'''
'''
		#Look for Pareto Distribution
		logM, logW = log_arrays(m,w)	
		
			
		SI = 0 #Start Index

		coeff = polyfit(logM[SI:],logW[SI:],1)
		print "power for line for $\\alpha = %.1f$: %f" % (0.5*(i+1),coeff[0])
		Tail = coeff[0]*logM[SI:] + coeff[1]
		
		plt.figure()
		plt.title("Pareto power %s" % f)
		plt.plot(logM[SI:],logW[SI:],".",markersize=0.8,label="$\\alpha = %s$" % 0.5*(i+1))
		plt.plot(logM[SI:],Tail,label="tail")
		plt.legend()
		plt.axis([-1,2,-5,0.5])
		plt.xlabel("$\log_{10}(m)$") ; plt.ylabel("$\log_{10}(w)$")

		#Look for Gibbs Distribution
		GSI = 0 ; GEI = 999	
		
		nlogM, nlogW = nlog_arrays(M[i][GSI:GEI],W[i][GSI:GEI]) ; m = exp(nlogM)
		
		Gcoeff = polyfit(m,nlogW,1)
		print "Gibbs slope for line %d : %f" % (i,Gcoeff[0])
		GTail = Gcoeff[0]*m + Gcoeff[1]
		
		plt.figure()
		plt.title("Gibbs slope %s" % f)
		plt.plot(m,nlogW,".",markersize = 0.8)
		plt.plot(m,GTail,label="Gtail")
		plt.legend()
		plt.xlabel("m") ; plt.ylabel("$\log(w)$")
		
	

	plt.figure()
	for i in range(5):
		logM, logW = log_arrays(M[i],W[i])
		plt.plot(logM,logW,".",markersize=0.8,label="%d" % i)
	plt.legend()
	plt.axis([-1,2,-5,0.5])
	plt.show()
'''
