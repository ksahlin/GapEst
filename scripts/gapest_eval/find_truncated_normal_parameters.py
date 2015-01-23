import argparse
from scipy.stats.distributions import norm
import math

def rekursion(mean,stddev,low_bound):
	"""
	Estimates the parameters of the normal distribution from truncated data

	"""
	mean_prev = mean
	while True:
		#print norm.pdf(low_bound,loc=mean, scale=stddev)
		low_bound_normalized = (low_bound - mean)/stddev
		lambda_low_bound = norm.pdf(low_bound_normalized) / (1-norm.cdf(low_bound_normalized))
		delta_low_bound = lambda_low_bound*(lambda_low_bound - low_bound_normalized)
		print lambda_low_bound, delta_low_bound
		mean = mean -  stddev*lambda_low_bound	
		var = stddev**2/(1-delta_low_bound ) 
		stddev = math.sqrt(var)
		print mean, stddev
		break
		if mean_prev - mean < 10:
			break 
		mean_prev = mean
	print mean, stddev
	return mean,stddev

def main(args):
	mean, stddev = rekursion(args.sample_mean, args.sample_stddev, args.low_bound)
if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('sample_mean', type=int, help='sample mean. ')
	parser.add_argument('sample_stddev', type=int, help='sample mean. ')
	parser.add_argument('low_bound', type=int, help='lower cutoff. ')

	args = parser.parse_args()
	main(args)