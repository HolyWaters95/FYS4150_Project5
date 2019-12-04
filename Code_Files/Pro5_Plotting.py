from Py_Functions import readarrays
from numpy import array, zeros, linspace, log, log10, sort, polyfit, polyval
import matplotlib.pyplot as plt
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
coeff = polyfit(m_intervals[:-1],w,10)
W = polyval(coeff,m_intervals[:-1])
plt.plot(m_intervals[:-1],W,'r')

fit_alpha, fit_loc, fit_beta=stats.gamma.fit(w)

plt.plot(
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

filenames = []
filestart = "../Results/Median"

Nvalues = ["500", "1000"]
Lvalues = ["0.250000", "0.500000", "0.900000"]
Avalues = ["0.500000", "1.000000", "1.500000", "2.000000"]
Gvalues = ["1.000000", "2.000000", "3.000000", "4.000000"]

#task a
filenames.append(filestart + "_N_" + Nvalues[0])

#task b
for 




