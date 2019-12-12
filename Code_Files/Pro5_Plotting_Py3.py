from Py_Functions_Py3 import readarrays, plot_median_c, nlog_arrays, plot_prob_distribution_c, log_arrays, f, non_zeros_array
import Py_Functions_Py3 as pf
from numpy import array, zeros, linspace, log, log10, exp, sort, polyfit, polyval
from matplotlib.pyplot import *
import matplotlib.image as pmg
import scipy.stats as stats



"""
task a)
"""
abc = input("Which task would you like to run)? a/b/c \n")
if abc == "a":
	N = range(1,501)
	Money = []
	M_R = readarrays("../Results_O/OE_a_Money_distributions_D_5_N_500_0.txt")[0]
	for i in range(0,len(M_R)):
		Money += list(M_R[i])

	w, m_intervals, patches = hist(Money,len(N),density=True)
	title("Histogram of 'amount of money' as function of $m$")
	xlabel("Amount of money")
	ylabel("$w_m$")
	savefig("../Plots/Plots_a_b_c/a_histogram.png", dpi=300)
	show()


"""
task b)
"""
if abc == "b":
	N = range(1,501)
	Money = []
	M_R = readarrays("../Results_O/OE_a_Money_distributions_D_5_N_500_0.txt")[0]
	for i in range(0,len(M_R)):
		Money += list(M_R[i])

	fig = figure()
	w, m_intervals, patches = hist(Money,len(N),density=True)
	close(fig)

	m = m_intervals[:-1]
	p, w = nlog_arrays(m, w)
	x, y = polyfit(exp(p[:-70]),w[:-70],1)

	plot(exp(p),w)
	plot(exp(p), x*exp(p)+y, label = "Slope = {0:f}%".format(x))
	title("Equilibrium, $w_m$, as function of $m$")
	xlabel("Wealth m (m_0 = 1)")
	ylabel("$log(w_m)$")
	legend()
	savefig("../Plots/Plots_a_b_c/b_logplot.png", dpi=300)
	show()






"""
task c
"""
if abc == "c":

	# Medianer
	D = 5
	filenames = []
	filestart = "../Results_O/OE_c_Median_D_%s" % D
	Nvalues = ["500"]
	Lvalues = ["0.250000", "0.500000", "0.900000"]
	for L in Lvalues:
		filenames.append(filestart + "_N_" + Nvalues[0] + "_L_" + L + ".txt")

	#for f in filenames:
	#	plot_median_c(f)


	# Distribution
	N = range(1,501)
	Money1 = []
	filenames = []
	filestart = "../Results_O/OE_c_Money_distributions_D_5_N_500"
	Lvalues = ["0.250000", "0.500000", "0.900000"]
	for L in Lvalues:
		filenames.append(filestart + "_L_" + L + "_0" + ".txt")
	for R in range(len(Lvalues)):
		Money1 = []
		M_R = readarrays("%s" % filenames[R])[0]
		for i in range(0,len(M_R)):
			Money1 += list(M_R[i])

		# Three plots of lambda and histograms
		figure()
		w, m_intervals, patches = hist(Money1,len(N),density=True)
		plot(m_intervals, f(Lvalues[R],m_intervals), label = "$\lambda = {0}$".format(Lvalues[R]))
		title("Money distributions for different $\lambda$")
		xlabel("Wealth m $(m_0 = 1)$")
		ylabel("$P(m)$")
		legend()
		savefig("../Plots/Plots_a_b_c/c_money_dist_hist" + "_L_" + Lvalues[R] + ".png", dpi=300)
		show()

	# All lambdas in the same plot
	for i in range(len(Lvalues)):
		plot(m_intervals, f(Lvalues[i],m_intervals), label = "$\lambda = {0}$".format(Lvalues[i]))
	title("Money distributions for different $\lambda$")
	xlabel("Wealth m $(m_0 = 1)$")
	ylabel("Percent of total participants $P(m)$")
	legend()
	savefig("../Plots/Plots_a_b_c/c_money_dist_L_123.png", dpi=300)
	show()


	# loglog plot, P(m) vs m
	for i in range(len(Lvalues)):
		plot(log(m_intervals), f(Lvalues[i],m_intervals), label = "$\lambda = {0}$".format(Lvalues[i]))
	title("Loglog plot of money distributions for different $\lambda$")
	xlabel("log(Wealth m $(m_0 = 1))$")
	ylabel("$log(P(m))$")
	yscale("log")
	legend()
	savefig("../Plots/Plots_a_b_c/c_money_dist_L_123_loglog.png", dpi=300)
	show()


	# plotting tail end
	Money2 = []
	filenames = []
	filestart = "../Results_O/OE_c_Money_distributions_D_5_N_500"
	Lvalues = ["0.250000"]
	for L in Lvalues:
		filenames.append(filestart + "_L_" + L + "_0" + ".txt")
	for R in range(len(Lvalues)):
		Money2 = []
		M_R = readarrays("%s" % filenames[R])[0]
		for i in range(0,len(M_R)):
			Money2 += list(M_R[i])

		fig = figure()
		w, m_intervals2, patches = hist(Money2,len(N),density=True)
		close(fig)
		m2 = m_intervals2[:-1]
		p2, w = nlog_arrays(m2, w)
		x2, y2 = polyfit(exp(p2[100:-50]),w[100:-50],1)

		plot(exp(p2[100:-50]),w[100:-50])
		plot(exp(p2[100:-50]), x2*exp(p2[100:-50])+y2, label = "Slope = {0:f}%".format(x2))

		#plot(p2, f(Lvalues[R],p2), label = "$\lambda = {0}$".format(Lvalues[R]))
		title("Equilibrium, $w_m$, as function of $m$")
		xlabel("Wealth m (m_0 = 1)")
		ylabel("$log(w_m)$")
		legend()
		savefig("../Plots/Plots_a_b_c/c_money_dist_L_123_tail.png", dpi=300)
		show()

		#plot(m_intervals, f(Lvalues[R],m_intervals), label = "$\lambda = {0}$".format(Lvalues[R]))







'''
coeff = polyfit(m_intervals[:-1],w,10)
W = polyval(coeff,m_intervals[:-1])
plt.plot(m_intervals[:-1],W,'r')

fit_alpha, fit_loc, fit_beta=stats.gamma.fit(w)

plt.plot()
'''



'''


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
"""
yn = input("Do you want to plot medians for task d)? y/n \n")
if yn == "y":
	filenames = []
	D = input("Exponent of D? \n")
	filestart = "../Results/Median_D_%s" % D
	for N in Nvalues:
		for L in ["0.000000","0.500000"]:
			for a in avalues:
				filenames.append(filestart + "_N_" + N + "_L_" + L + "_a_" + a + ".txt")

	for f in filenames:
		plot_median_d(f,save=False)

yn = input("Do you want to plot probability distributions for d)? y/n \n")
if yn == "y":
	filenames = []
	D = input("Exponent of D? \n")
	filestart = "../Results/Money_distributions_D_%s" % D
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

		pf.plot_loglogW(M[i:i+4],W[i:i+4],filenames[i:i+4],save=True)
		pf.Pareto_dist_d(M[i:i+4],W[i:i+4],filenames[i:i+4],save=True)
		pf.Gibbs_dist_d(M[i:i+4],W[i:i+4],filenames[i:i+4],save=True)
		i += 4

yn = input("Do you want to plot medians for task e)? y/n \n")
if yn == "y":
	filenames = []
	D = input("Exponent of D? \n")
	filestart = "../Results/Median_D_%s_N_1000" % D
	for L in ["0.000000","0.500000"]:
		for a in ["1.000000","2.000000"]:
			for g in gvalues:
				filenames.append(filestart + "_L_" + L + "_a_" + a + "_g_" + g + ".txt")
	filenames.append("../Results/Median_D_8_N_1000_L_0.500000_a_2.000000_g_1.000000.txt")
	#for f in filenames:
	#	pf.plot_median_e(f,save=True)
	pf.plot_median_e(filenames[-1],save=True)

w,m = plot_prob_distribution_d("../Results/Money_distributions_D_7_N_1000_L_0.000000_a_2.000000_g_2.000000.txt")

logm, logw = log_arrays(m[:-1],w)
plot(logm,logw)
show()
"""
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
