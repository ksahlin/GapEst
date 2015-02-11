import argparse
from scipy.stats.distributions import norm
#import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from scipy.optimize import fmin
import random
from assemblymodule import lib_est


def estimate_params_for_normal(x, low_bound , mu_initial, sigma_initial):
	"""
		Takes a vector x of truncated data with a known lower
		truncation bound and estimates the parameters of the 
		fit of an untruncated normal distribution.
		code from Chris Fonnesbeck's Python data analysis tutorial on Sense
		https://sense.io/prometheus2305/data-analysis-in-python/files/Statistical%20Data%20Modeling.py
	"""


	# normalize vector
	mu_initial = float(mu_initial)
	sigma_initial = float(sigma_initial)
	#x = np.random.normal(size=10000,loc=2000,scale= 2000)

	x = map(lambda y: (y-mu_initial )/sigma_initial ,x)
	a =  (low_bound - mu_initial)/sigma_initial # normalize lower bound
	

	#_ = plt.hist(x, bins=100)
	#plt.show()
	#plt.close()

	# We can construct a log likelihood for this function using the conditional
	# form	
	trunc_norm = lambda theta, a, x: -(np.log(norm.pdf(x, theta[0], theta[1])) - 
	                                      np.log(1 - norm.cdf(a, theta[0], theta[1]))).sum()

	# For this example, we will use another optimization algorithm, the
	# **Nelder-Mead simplex algorithm**. It has a couple of advantages: 
	# 
	# - it does not require derivatives
	# - it can optimize (minimize) a vector of parameters
	# 
	# SciPy implements this algorithm in its `fmin` function:

	# we have normalized data, given that the loer truncation point a
	# is pretty far out in the tail - the standard normal parameters are
	# a first good guess, i.e. 0,1
	initial_guess = np.array([0,1]) 
	sol = fmin(trunc_norm, initial_guess , args=(a, x))
	print sol
	mean_normalized,stddev_normalized = sol[0],sol[1]
	mean_est =( 1 + mean_normalized ) * mu_initial
	stddev_est = stddev_normalized * sigma_initial
	print mean_est,stddev_est
	return mean_est,stddev_est


# def rekursion(mean,stddev,low_bound):
# 	"""
# 	Estimates the parameters of the normal distribution from truncated data

# 	"""
# 	# work with normalized values
# 	low_bound_normalized = (low_bound - mean)/stddev
# 	print low_bound_normalized
# 	mean_normalized =  0
# 	stddev_normalized = 1

# 	# find converged parameters
# 	mean_prev = mean_normalized
# 	i=0
# 	while i< 10:
# 		lambda_low_bound = norm.pdf(low_bound_normalized) / (1-norm.cdf(low_bound_normalized))
# 		delta_low_bound = lambda_low_bound*(lambda_low_bound - low_bound_normalized)
# 		#print lambda_low_bound, delta_low_bound
		
# 		mean_normalized = mean_normalized -  stddev_normalized*lambda_low_bound	
# 		var = stddev_normalized**2/(1-delta_low_bound ) 
# 		stddev_normalized = math.sqrt(var)
# 		low_bound_normalized = 000000000000#
# 		print mean_normalized, stddev_normalized
# 		# break
# 		if mean_prev - mean_normalized < 0.001:
# 			break 
# 		mean_prev = mean_normalized
# 		i+=1

# 	# transform back
# 	est_mean = ( 1 - mean_normalized)*mean
# 	est_stddev = (stddev_normalized)*stddev
# 	print est_mean,est_stddev 
# 	return est_mean,est_stddev
def get_bam_observations(bamfile_path):
	#read_length = 101
	low_bound =  202
	#bam_object = CreateGraph_updated.BamParser(bamfile_path)
	i_sizes = []
	i = 0
	for isize in lib_est.proper_read_isize_iter(bamfile_path,101,12000): #bam_object.aligned_reads('bwa'):
	    i_sizes.append(isize)
	    # if read.is_read1 and CreateGraph_updated.is_proper_aligned_unique_innie(read):
	    # 	if abs(read.tlen) > low_bound: 
	    #     	i_sizes.append(abs(read.tlen))
	    #     	i+=1
	    #if i > 100000:
	    #	break 

	print 'Nr proper reads:', len(i_sizes)
	filtered_observations = GapCalculator.remove_misalignments(i_sizes,10)
	n_isize = float(len(filtered_observations))
	mean_isize = sum(filtered_observations) / n_isize
	std_dev_isize = (sum(list(map((lambda x: x ** 2 - 2 * x * mean_isize + mean_isize ** 2), filtered_observations))) / (n_isize - 1)) ** 0.5
	# low_bound = 2.0*read_length
	return filtered_observations, mean_isize, std_dev_isize,low_bound

def get_bp_observations(bp_file):
	x = []

	for i, line in enumerate(open(bp_file,'r')):
		scf, pos, n_obs, mean, stddev = line.strip().split()
		pos, n_obs,mean,stddev = int(pos), float(n_obs), float(mean), float(stddev)
		if float(n_obs) < 2 or float(mean) < 0:
			continue
		else:
			x.append(mean)
		if i % 4000001 == 0:
			break

	n = float(len(x))
	mu = sum(x)/n
	sigma = (sum(list(map((lambda t: t ** 2 - 2 * t * mu + mu ** 2), x))) / (n - 1)) ** 0.5
	return x, mu, sigma, 0

def main(args):
	if args.bamfile:
		x, mu, sigma,low_bound = get_bam_observations(args.bamfile)
	elif args.bp_file:
		x, mu, sigma,low_bound = get_bp_observations(args.bp_file)

	# x_reduced = random.sample(x, 10000)
	# print min(x_reduced), max(x_reduced)

	mu_est, sigma_est = estimate_params_for_normal(x, low_bound, mu, sigma )

	#plt fit
	_ = plt.hist(x, bins=200, normed=True, alpha=0.6, color='g')
	xmin, xmax = plt.xlim()
	#mu, sigma = norm.fit(x)
	#print mu,sigma
	y = np.linspace(xmin, xmax, 100)
	plt.plot(y,mlab.normpdf(y,mu_est, sigma_est),'r')
	title = "Fit results: mu = %.2f,  std = %.2f" % (mu_est, sigma_est)
	plt.title(title)
	plt.show()
	#mean, stddev = rekursion(args.sample_mean, args.sample_stddev, args.low_bound)
if __name__ == '__main__':
	from src import CreateGraph_updated, GapCalculator
	parser = argparse.ArgumentParser()
	parser.add_argument('--bam_file', dest='bamfile',type=str, default=False, help='bamfile. ')
	parser.add_argument('--bp_file', dest='bp_file', type=str, default = False, help='Base pair stats file generated by GetDistr-assembly module. ')

	args = parser.parse_args()
	main(args)
