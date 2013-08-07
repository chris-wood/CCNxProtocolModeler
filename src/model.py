# File: ccn_model.py
# Author: Christopher A. Wood

import sys
from math import *
import matplotlib.pyplot as plt
import time

# References:
# [1] - Modeling data transfer in content-centric networking (extended version)

# Model parameters used to specify constraints on the topology/scenario being modeled
rts = 0
N = 0
M = 0
K = 0
c = 0
x = 0
alpha = 0
sigma = 0
sigma_k = []
delta = 0
delta_k = []
lmbda = 0
lmbda_k = []

### MODEL EQUATIONS ###

def calc_g():
	''' Calculate g according to equation (2) in [1]

		Status: Correct.
	'''
	global lmbda, alpha, sigma, M, K, c

	m = float(M) / float(K)
	g = lmbda * c * pow(sigma, alpha) * pow(m, (alpha - 1.0)) * pow(gamma(1.0 - (1.0 - alpha)), alpha)
	return 1.0 / g

def calc_q(k):
	global c, alpha
	return c / (k ** alpha)

def calc_p1(k):
	''' Calculate p_k(1) according to first part of proof of proposition 6.2 in [1]

	 	p_k(1) = e^( (-lambda / m) * q_k * g * x^alpha )

		Status: Correct.
	'''
	global lmbda, alpha, M, K, c, x

	m = float(M) / float(K)
	q_k = calc_q(k)
	g = calc_g()
	p_k1 = exp( ( (lmbda * -1.0) / m ) * q_k * pow(x, alpha) * g )
	return p_k1

def calc_p_filter(k, i, filtering = True):
	''' Calculate p_filt(k,i) according to equation (28) in [1]

		Status: Correct.
	'''
	global lmbda, alpha, sigma, sigma_k, delta, delta_k, M, K, c, x

	if filtering: # for k \in {1,..,N}
		bki = calc_bk(k, i, filtering)
		bki_sigma = pow(bki, sigma_k[k - 1])
		num = bki * ( 1.0 - bki_sigma )
		denom = (1.0 - bki) * sigma_k[k - 1]
		frac = num / denom
		return (1.0 - frac)
	else:
		raise Exception("Non-filtering model not implemented.")

def calc_bk(k, i, filtering = True):
	''' Calculate b_k(i) according to equation 28 in [1]

		Status: Match.
	'''
	global lmbda, alpha, sigma, sigma_k, delta, delta_k, M, K, c, x

	if filtering:
		m = float(M / K)
		if (i == 1):
			return exp( (delta_k[k - 1] * -1) * sigma_k[k - 1] / m )
		else:
			mu_k = calc_mu(k, i - 1, filtering)
			return exp( (delta_k[k - 1] * -1) * 2.0 * mu_k / m)
	else:
		raise Exception("Non-filtering model not implemented.")

def calc_mu(k, i, filtering = True):
	''' Calculate mu_k(i) according to equation 28 in [1]

		Status: Match.
	'''	
	global lmbda, lmbda_k, alpha, sigma, sigma_k, delta, delta_k

	if filtering:
		if (i == 1):
			return (lmbda_k[k - 1] * calc_p1(k) * (1.0 - calc_p_filter(k, 1, filtering)))
		else:
			return (2.0 * calc_mu(k, i - 1, filtering) * calc_p(k, i, filtering) * (1.0 - calc_p_filter(k, i, filtering)))
	else:
		raise Exception("Non-filtering model not implemented.")

def calc_p(k, i, filtering = True): # assume filtering (interest aggregation) is enabled
	''' Calculate p_k(i) according to proposition 6.5 in [1].

		p_k(1) = e^([-\lambda / m)] * q_k * g * x^alpha)

		Status: 
	'''
	global lmbda, alpha, sigma, M, K, c, x
	
	# M = number of different content items, K = # of content classes, so there's m things in each class k
	m = float(M) / float(K)

	if filtering:
		pk1 = calc_p1(k) # p_k(1)^f \equiv p_k(1) (no previous filtering has occurred)

		# Sanity check
		if (pk1 > 1.0):
			raise Exception("Pk(1) can't be larger than 1.0")

		# Calculate the exponent
		exponent = float(1.0) 
		for l in range(1, i): #[1, i - 1]
			p_k_l = calc_p(k, l, filtering)
			p_filter_k_l = calc_p_filter(k, l)

			# More sanity checks
			if (p_k_l > 1.0):
				raise Exception("ERROR: miss probability invalid: " + str(p_k_l))
			if (p_filter_k_l > 1.0):
				raise Exception("ERROR: filter miss invalid: " + str(p_filter_k_l))

			# Accumulate the exponent
			# print("Round: " + str(l))
			# print(p_k_l)
			# print(p_filter_k_l)
			# print("pre exp: " + str(exponent))

			# FROM TECHNICAL REPORT
			# exponent = exponent * (p_k_l * (1.0 - p_filter_k_l))

			# FROM PAPER
			exponent = exponent * p_k_l 


			# print("post exp: " + str(exponent))

		# p = p_k(1) ^ exponent
		# print("base, exp")
		# print(pk1)
		# print(exponent)
		# print(pow(pk1, exponent))
		# print("done")
		# time.sleep(1)
		return pow(pk1, exponent)
	else:
		raise Exception("Non-filtering model not implemented.")

def calc_rt(i):
	global rts
	return 2.0 * rts * i

# the model assumes that all links have the same round trip delay, which is fine
# this is the virtual RTT for a single chunk, which is ultimately part of a piece of content somehow, depending on fragmentation

# Status: Validated.
def calc_VRTT(k):
	global N
	global rts
	s = 0.0

	# Compute the weighted sum
	for i in range(1, N + 1): #[1, N]
		pki = calc_p(k, i)
		# print("k,i,pki = " + str(k) + ","+ str(i) + "," + str(pki))
		if (pki > 1.0):
			raise Exception("Probabilities (pki) can't be larger than one")
		rti = calc_rt(i) # $R_i$
		prod = float(1.0)
		for j in range(1, i): #[1, i - 1]
			pkj = calc_p(k, j)
			if (pkj > 1.0):
				raise Exception("Probabilities (pkj) can't be larger than one")
			prod = prod * pkj
		s = s + (rti * (1 - pki) * prod)

	# Safety check
	if (s < 0):
		print("k = " + str(k))
		raise Exception("Can't have negative VRTT")

	return s

def calcAverageThroughput(W, k):
	return W / calc_VRTT(k)

def calcTime(sigma, W, k):
	X = calcAverageThroughput(W, k)
	return sigma / (W * X)

def main():
	# Read parameters in from a file and calculate accordingly
	global N 
	global M 
	global rts 
	global K 
	global alpha 
	global c 
	global x 
	global sigma 
	global sigma_k 
	global lmbda 
	global lmbda_k 
	global delta 
	global delta_k 

	print("Populating with some test data...")

	# Taken from the paper... let's see if this works.
	M = 20000
	K = 400
	alpha = float(1.2)
	lmbda = float(40)
	W = 1.0
	x = float(200000) # 2GBs
	sigma = float(690) # 10Kb chunks on average... but ideally, we'd like to model this using some random distribution
	c = 1.0 
	N = 3 # levels of the binary tree
	rts = float(2) # 2ms round trip time between nodes
	delta = float(1) 

	# Generate the data sets for each specific class of content items...
	sigma_k = []
	delta_k = [] # technically, delta_k = min(delta, VRTT_k), but we need delta_k to compute VRTT_k... cycle.
	lmbda_k = []
	for i in range(K):
		sigma_k.append(sigma)
		delta_k.append(delta)
		# lmbda_k.append(lmbda)
		lmbda_k.append(lmbda * calc_q(i + 1))

	# Do the calculations
	n1 = []
	mu1 = []
	vrtt1 = []
	n2 = []
	mu2 = []
	vrtt2 = []
	n3 = []
	mu3 = []
	vrtt3 = []
	vrtt = []
	for k in range(1, K + 1):
		vrtt.append(calc_p(k, 2))
		# mu1.append(calc_mu(k, 1, True))
		# mu2.append(calc_mu(k, 2, True))
		# n1.append(calc_p(k, 1))
		# n2.append(calc_p(k, 2))
		# n3.append(calc_p(k, 3))

	# Show the plot
	x = range(K)
	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	# ax1.scatter(x, n1, s=10, c='b', marker="s")
	# ax1.scatter(x, vrtt, s=10, c='r', marker="s")
	# ax1.scatter(x, mu1, s=10, c='g', marker="s")
	ax1.scatter(x, vrtt, s=10, c='b', marker="s")
	plt.show()

	# for v in lmbda_k:
	# 	print(v)

if __name__ == "__main__":
	main()