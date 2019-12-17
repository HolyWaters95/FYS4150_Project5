from Py_Functions_Py3 import readarrays, plot_median_c, nlog_arrays, plot_prob_distribution_c, log_arrays, P, non_zeros_array
import Py_Functions_Py3 as pf
from numpy import array, zeros, linspace, log, log10, exp, sort, polyfit, polyval
from matplotlib.pyplot import *
import matplotlib.image as pmg
import scipy.stats as stats



"""
task a)
"""
abc = input("Which task would you like to run? a/b/c \n")
if (abc == "a"):
	N = range(1,501)
	Money = []
	M_R = readarrays("../Results_O/OE_a_Money_distributions_D_5_N_500_0.txt")[0]
	for i in range(0,len(M_R)):
		Money += list(M_R[i])

	w, m_intervals, patches = hist(Money,len(N)*2,density=True, label="Numerical")
	plot(m_intervals, exp(-m_intervals), label = "Analytical", linewidth = 0.9)
	title("Wealth distribution with $N$ = 500")
	xlabel("Wealth m ")
	ylabel("$w_m$")
	legend()
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
	x, y = polyfit(exp(p[:-60]),w[:-60],1)

	plot(exp(p),w)
	plot(exp(p), x*exp(p)+y, label = "Slope = {0:.4}".format(x))
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
	"""
	# Medianer
	D = 5
	filenames = []
	filestart = "../Results_O/OE_c_Median_D_%s" % D
	Nvalues = ["500"]
	Lvalues = ["0.250000", "0.500000", "0.900000"]
	for L in Lvalues:
		filenames.append(filestart + "_N_" + Nvalues[0] + "_L_" + L + ".txt")

	for f in filenames:
		plot_median_c(f)


	# Distribution
	N = range(1,501)
	Money = []
	filenames = []
	filestart = "../Results_O/OE_c_Money_distributions_D_5_N_500"
	Lvalues = ["0.250000", "0.500000", "0.900000"]
	for L in Lvalues:
		filenames.append(filestart + "_L_" + L + "_0" + ".txt")
	for R in range(len(Lvalues)):
		Money = []
		M_R = readarrays("%s" % filenames[R])[0]
		for i in range(0,len(M_R)):
			Money += list(M_R[i])

		# Three seperate plots of lambda superimposed on histograms
		figure()
		w, m_intervals, patches = hist(Money,len(N),density=True, label="Numerical")
		plot(m_intervals, P(Lvalues[R],m_intervals), label = "Analytical")
		title("Distribution of wealth for $\lambda = ${0:.2f}".format(float(Lvalues[R])))
		xlabel("Wealth m $(m_0 = 1)$")
		ylabel("$P(m)$")
		legend()
		savefig("../Plots/Plots_a_b_c/c_money_dist_hist" + "_L_" + Lvalues[R] + ".png", dpi=300)
		show()
		close()


	# All lambdas in the same plot
	# Analytical lambda = 0
	N = range(1,501)
	Money1 = []
	M_R = readarrays("../Results_O/OE_a_Money_distributions_D_5_N_500_0.txt")[0]
	for i in range(0,len(M_R)):
		Money1 += list(M_R[i])
	fig = figure()
	un_used, m_intervals0, patches = hist(Money1,len(N),density=True, label="Numerical")
	close(fig)
	plot(m_intervals0[14:-417], exp(-m_intervals0[14:-417]), label = "Analytical $\lambda = 0$")
	# Analytical Lambda = 0.25, 0.5, 0.9
	for i in range(len(Lvalues)):
		plot(m_intervals, P(Lvalues[i],m_intervals), label = "Analytical $\lambda = ${0:.2f}".format(float(Lvalues[i])))
	title("Distribution of wealth for different $\lambda$")
	xlabel("Wealth m $(m_0 = 1)$")
	ylabel("Percent of total participants $P(m)$")
	legend()
	savefig("../Plots/Plots_a_b_c/c_money_dist_L_123.png", dpi=300)
	show()


	# loglog plot of all lambda
	# Lambda = 0
	plot(log(m_intervals), exp(-m_intervals), label = "Analytical $\lambda = 0$")
	# Lambda = 0.25, 0.5, 0.9
	for i in range(len(Lvalues)):
		plot(log(m_intervals), P(Lvalues[i],m_intervals), label = "Analytical $\lambda = ${0:.2f}".format(float(Lvalues[i])))
	title("Loglog plot of Distribution of wealth for different $\lambda$")
	xlabel("log(Wealth m $(m_0 = 1))$")
	ylabel("$log(P(m))$")
	yscale("log")
	legend()
	savefig("../Plots/Plots_a_b_c/c_money_dist_L_123_loglog.png", dpi=300)
	show()
	"""

	# plotting tail end for each lambda
	N = range(1,501)
	filenames = []
	filestart = "../Results_O/OE_c_Money_distributions_D_5_N_500"
	Lvalues = ["0.250000", "0.500000", "0.900000"]
	for L in Lvalues:
		filenames.append(filestart + "_L_" + L + "_0" + ".txt")
	for R in range(len(Lvalues)):
		Money2 = []
		M_R = readarrays("%s" % filenames[R])[0]
		for i in range(0,len(M_R)):
			Money2 += list(M_R[i])
		fig = figure()
		ww, m_intervals2, patches = hist(Money2,len(N),density=True)
		close(fig)
		m2 = m_intervals2[:-1]
		p2, www = nlog_arrays(m2, ww)

		start = (150, 200, 300)
		stop = (-40, -40, -10)
		# tail for lambda = 0.25 with fitting
		if R == 0:
			x1, y1 = polyfit(exp(p2[start[0]:stop[0]]),www[start[0]:stop[0]],1)
			plot(exp(p2[start[0]:stop[0]]),www[start[0]:stop[0]])
			plot(exp(p2[start[0]:stop[0]]), x1*exp(p2[start[0]:stop[0]])+y1, label = "Slope = {0:.3f}".format(float(x1)))
			title("Tail end of distribution of wealth $\lambda$ = {0:.2f}".format(float(Lvalues[0])))
			xlabel("Wealth m (m_0 = 1)")
			ylabel("$log(w_m)$")
			legend()
			axis(0,-10,1,6)
			savefig("../Plots/Plots_a_b_c/c_money_dist_L_" + Lvalues[0] + "_tail.png", dpi=300)
			show()
		# tail for lambda = 0.50 with fitting
		if R == 1:
			x2, y2 = polyfit(exp(p2[start[1]:stop[1]]),www[start[1]:stop[1]],1)
			plot(exp(p2[start[1]:stop[1]]),www[start[1]:stop[1]])
			plot(exp(p2[start[1]:stop[1]]), x2*exp(p2[start[1]:stop[1]])+y2, label = "Slope = {0:.3f}".format(float(x2)))
			title("Tail end of distribution of wealth $\lambda$ = {0:.2f}".format(float(Lvalues[1])))
			xlabel("Wealth m (m_0 = 1)")
			ylabel("$log(w_m)$")
			legend()
			axis(0,-10,1,6)
			savefig("../Plots/Plots_a_b_c/c_money_dist_L_" + Lvalues[1] + "_tail.png", dpi=300)
			show()
		# tail for lambda = 0.90 with fitting
		if R == 2:
			x3, y3 = polyfit(exp(p2[start[2]:stop[2]]),www[start[2]:stop[2]],1)
			plot(exp(p2[start[2]:stop[2]]),www[start[2]:stop[2]])
			plot(exp(p2[start[2]:stop[2]]), x3*exp(p2[start[2]:stop[2]])+y3, label = "Slope = {0:.3f}".format(float(x3)))
			title("Tail end of distribution of wealth $\lambda$ = {0:.2f}".format(float(Lvalues[2])))
			xlabel("Wealth m (m_0 = 1)")
			ylabel("$log(w_m)$")
			legend()
			axis(0,-10,1,6)
			savefig("../Plots/Plots_a_b_c/c_money_dist_L_" + Lvalues[2] + "_tail.png", dpi=300)
			show()
