'''
Created on Sep 18, 2013

@author: ksahlin
'''

import argparse
import os
from mathstats.normaldist.normal import MaxObsDistr
from scipy.stats import ks_2samp,norm
import random
#import re
import math

import pysam
#import math
from itertools import ifilter
#import model
import bisect

from collections import deque

#import matplotlib.pyplot as plt
from statsmodels.distributions.empirical_distribution import ECDF

EMPIRICAL_BINS = 200
SAMPLE_SIZE = 50000  # for estimating true full read pair distribution

def is_proper_aligned_unique_innie(read):
    return (read.is_reverse and not read.mate_is_reverse and read.is_read1 and read.tlen < 0 and read.rname == read.mrnm) or \
                (not read.is_reverse and read.mate_is_reverse and read.is_read1 and read.tlen > 0 and read.rname == read.mrnm ) \
                and not read.mate_is_unmapped and read.mapq > 10 and not read.is_secondary

def ReadInContigseqs(contigfile,contig_filter_length):
    cont_dict = {}
    k = 0
    temp = ''
    accession = ''
    for line in contigfile:
        if line[0] == '>' and k == 0:
            accession = line[1:].strip().split()[0]
            cont_dict[accession] = ''
            k += 1
        elif line[0] == '>':
            cont_dict[accession] = temp
            temp = ''
            accession = line[1:].strip().split()[0]
        else:
            temp += line.strip()
    cont_dict[accession] = temp

    # for contig in cont_dict.keys():
    #     print 'Initial length:', len(cont_dict[contig])
    if contig_filter_length:
        singled_out = 0
        for contig in cont_dict.keys():
            if len(cont_dict[contig]) < contig_filter_length:
                del cont_dict[contig]
                singled_out += 1
    return(cont_dict)


class Parameters(object):
	"""docstring for Parameters"""
	def __init__(self):
		super(Parameters, self).__init__()
		self.mean = None
		self.stddev = None
		self.d = None
		self.pval = None
		self.scaf_lengths = {}
		self.nobs = None
		self.true_distr = None
		self.corrected_pval = None

	def get_pval_threshold(self):
		mean_sd = t.ppf( 0.975, self.nobs - 1 ) * self.stddev / math.sqrt( self.nobs )
		sd_sd = math.sqrt( (self.nobs - 1) * self.stddev**2 / chi2.ppf( 0.025, self.nobs - 1 ) )

		mean_conservative = self.mean + mean_sd
		sd_conservative = sd_sd

		z = norm.ppf( self.pval / 2.0, mean_conservative, sd_conservative )
		return norm.cdf( z, self.mean, self.stddev )

	def sample_distribution(self,bamfile,outfile):
		isize_list = []
		#i = 0
		bam_filtered = ifilter(lambda r: is_proper_aligned_unique_innie(r), bamfile)
		#while i <= sample_size:
		for sample_nr,read in enumerate(bam_filtered):
	   		## add do insert size distribution calculation if proper pair
			if is_proper_aligned_unique_innie(read):
				isize_list.append(abs(read.tlen))
				#sample_nr+=1
			if sample_nr > SAMPLE_SIZE:
				break
		print >> outfile, '#Insert size sample size:', sample_nr
		bamfile.reset()

		n_isize = float(len(isize_list))
		mean_isize = sum(isize_list)/n_isize
		std_dev_isize =  (sum(list(map((lambda x: x ** 2 - 2 * x * mean_isize + mean_isize ** 2), isize_list))) / (n_isize - 1)) ** 0.5
		print >> outfile,'#Mean before filtering :', mean_isize
		print >> outfile,'#Stddev before filtering: ', std_dev_isize
		extreme_obs_occur = True
		while extreme_obs_occur:
			extreme_obs_occur, filtered_list = AdjustInsertsizeDist(mean_isize, std_dev_isize, isize_list)
			n_isize = float(len(filtered_list))
			mean_isize = sum(filtered_list) / n_isize
			std_dev_isize = (sum(list(map((lambda x: x ** 2 - 2 * x * mean_isize + mean_isize ** 2), filtered_list))) / (n_isize - 1)) ** 0.5
			isize_list = filtered_list

		print >> outfile,'#Mean converged:', mean_isize
		print >> outfile,'#Std_est converged: ', std_dev_isize
		print >> outfile,'{0}\t{1}'.format( mean_isize, std_dev_isize)
		self.nobs = n_isize
		self.mean = mean_isize
		self.stddev = std_dev_isize 
		self.full_ECDF = ECDF(isize_list)
		self.adjustedECDF_no_gap = None
		self.get_true_normal_distribution(random.sample(isize_list, min(10000,sample_nr)),outfile)


	def get_true_normal_distribution(self,sample,outfile):
		read_len = 100
		softclipps = 0
		#self.temp_cdf = {}
		#self.temp_cdf[norm.cdf(2*(read_len-softclipps), self.mean, self.stddev) * 1] = 2*(read_len-softclipps)  # weight one for one placement
		cdf_list =[norm.cdf(2*(read_len-softclipps), self.mean, self.stddev) * 1]
		for x in range(2*(read_len-softclipps) +1,int(self.mean + 6*self.stddev)):
			increment_area = norm.pdf(x, self.mean, self.stddev) * (x-(2*(read_len-softclipps)-1))
			#self.temp_cdf[sum(cdf_list) + increment_area] = x
			cdf_list.append( cdf_list[-1] + increment_area)
		tot_cdf = cdf_list[-1]
		cdf_list_normalized = map(lambda x: x /float(tot_cdf),cdf_list)

		# Now create a weighted sample
		self.true_distr = []
		for i in range(1000):
			obs = random.uniform(0, 1)
			pos = bisect.bisect(cdf_list_normalized, obs) - 1
			#print obs, pos
			self.true_distr.append(pos + 2*(read_len-softclipps))


		n = len(self.true_distr)
		self.adjusted_mean = sum(self.true_distr)/float(len(self.true_distr))
		self.adjusted_stddev = (sum(list(map((lambda x: x ** 2 - 2 * x * self.adjusted_mean + self.adjusted_mean ** 2), self.true_distr))) / (n - 1)) ** 0.5

		print >> outfile,'#Corrected mean:{0}, corrected stddev:{1}'.format(self.adjusted_mean, self.adjusted_stddev)
		print >> outfile,'{0}\t{1}'.format(self.adjusted_mean, self.adjusted_stddev)




	def get_weight(self,x,gap_coordinates,r,s):
		if not gap_coordinates:
			return x - (2*(r-s)-1)
		total_restriction_positions_left = 0
		total_restriction_positions_right = 0

		for start,stop in gap_coordinates:
			if  x < start:
				continue #total_restriction_positions_right += 0

			elif  (r-s)-1 <= start  <= x:
				#print '1'
				total_restriction_positions_right += min(stop + (r-s)-1, x) - start  + (r-s)-1
			
			elif 0 <= start <= (r-s)-1:
				#print '2'
				total_restriction_positions_right += min(stop + (r-s)-1, x) 

			elif stop < -x:
				total_restriction_positions_left += (r-s)

			elif -x <= stop < -( (r-s)-1):
				#print '3'
				total_restriction_positions_left += stop - max(start - (r-s)-1 ,-x)  + (r-s)-1
				
			elif -( (r-s)-1) <= stop <= 0:
				#print '4'
				total_restriction_positions_left += -(max(start,-x))

			elif start <0 and stop > 0:
				#print '5'
				total_restriction_positions_right +=  stop + (r-s)
				total_restriction_positions_left += - start + (r-s)

		#print 'tot restrict:', total_restriction_positions
		total_restriction_positions_right = max(total_restriction_positions_right,s)
		total_restriction_positions_left = max(total_restriction_positions_left,s)
		weight = x - total_restriction_positions_right - total_restriction_positions_left

		return max( 0 , weight)



			# elif stop <= 0 and -stop < x < -start:
			# 	pass
			# elif start <0 and stop > 0:
			# 	pass
		
		return 

	def get_correct_ECDF(self,outfile, gap_coordinates):
		if not gap_coordinates and self.adjustedECDF_no_gap:
			return self.adjustedECDF_no_gap


		read_len = 100
		softclipps = 80


		x_min = max(2*(read_len-softclipps) , int(self.mean - 5*self.stddev) )
		x_max = int(self.mean + 5*self.stddev)
		stepsize =  (x_max - x_min) / EMPIRICAL_BINS
		cdf_list = [ self.full_ECDF( x_min) * self.get_weight(x_min, gap_coordinates, read_len, softclipps)  ] #[ self.full_ECDF( 2*(read_len-softclipps)) * self.get_weight(2*(read_len-softclipps), gap_coordinates, read_len, softclipps) ]


		# create weigted (true) distribution

		for x in range( x_min + stepsize , x_max, stepsize):
			increment_area = self.get_weight(x,gap_coordinates, read_len, softclipps) * (self.full_ECDF(x) - self.full_ECDF(x-stepsize))
			cdf_list.append( cdf_list[-1] + increment_area)
		# cdf_list = [self.get_weight(x,gap_coordinates, read_len, softclipps) * (self.full_ECDF(x) - self.full_ECDF(x-stepsize)) ]

		#print 'stepsize:', stepsize
		#print 'BINS:',len(cdf_list)
		tot_cdf = cdf_list[-1]
		cdf_list_normalized = map(lambda x: x /float(tot_cdf),cdf_list)

		# Now create a weighted sample

		self.true_distr = [ (bisect.bisect(cdf_list_normalized, random.uniform(0, 1) ) - 1)*stepsize + x_min for i in range(1000) ]

		# initialization of no gap true distribution
		if not gap_coordinates:
			print 'getting initial gap free distr.'
			self.adjustedECDF_no_gap = self.true_distr

		n = len(self.true_distr)
		self.adjusted_mean = sum(self.true_distr)/float(len(self.true_distr))
		self.adjusted_stddev = (sum(list(map((lambda x: x ** 2 - 2 * x * self.adjusted_mean + self.adjusted_mean ** 2), self.true_distr))) / (n - 1)) ** 0.5

		#print 'Corrected mean:{0}, corrected stddev:{1}, gap_coordinates: {2}'.format(self.adjusted_mean, self.adjusted_stddev, gap_coordinates)

		#print >> outfile,'Corrected mean:{0}, corrected stddev:{1}, gap_coordinates: {2}'.format(self.adjusted_mean, self.adjusted_stddev, gap_coordinates)
		return self.true_distr

class ReadContainer(object):
	"""docstring for ReadContainer"""
	def __init__(self, position):
		super(ReadContainer, self).__init__()
		self.position = position
		self.reads = []
	def __str__(self):
		"Prints reads on fasta format"
		pass
	def add_read(self,read):
		self.reads.append(read)

	def print_bam(self):
		"Print reads on bam format "


	def write_pval_to_file(self,outfile,ref_name):
		print >> outfile, '{0}\t{1}\t{2}\t{3}\t{4}\t{5}'.format(ref_name, self.position, self.pval,len(self.isize_list), self.mean_isize, self.std_dev_isize)

	def calc_observed_insert(self):
		#print self.reads
		#print 'list:',map(lambda x: x.tlen, self.reads )
		self.mean_isize = 0
		self.std_dev_isize = 0
		self.isize_list = map(lambda x: x.tlen, self.reads )
		#print len( self.isize_list)
		if len( self.isize_list) <= 5:
			self.mean_isize = -1
			return

		n_isize = float(len(self.isize_list))
		self.mean_isize = sum(self.isize_list)/n_isize
		self.std_dev_isize =  (sum(list(map((lambda x: x ** 2 - 2 * x * self.mean_isize + self.mean_isize ** 2), self.isize_list))) / (n_isize - 1)) ** 0.5



	def calc_ks_test(self,true_distribution):
		if len(self.isize_list) >= 5:
			KS_statistic, self.pval = ks_2samp(self.isize_list, true_distribution)
			return KS_statistic, self.pval 
		else:
			self.pval = -1
			return -1, -1


class BreakPointContainer(object):
	"""docstring for BreakPointContainer"""
	def __init__(self,param):
		super(BreakPointContainer, self).__init__()
		self.clusters = {}
		self.index = 1
		self.clusterinfo = {}
		self.param = param

	def add_bp_to_cluster(self, scf, pos, p_val, nr_obs, mean_obs, sv_type_observed, window_size):
		new_cluster = 1
		for (sv_type, ref, i), cluster in self.clusters.iteritems():
			if scf != ref:
				continue
			min_pos = cluster['min']
			max_pos = cluster['max']

			if sv_type_observed != sv_type:
				continue

			if pos <= max_pos + window_size and pos >= min_pos - window_size:
				new_cluster = 0
				cluster['pos'].append(pos)
				if pos < min_pos:
					 cluster['min'] = pos
				if pos > max_pos:
					 cluster['max'] = pos

				self.clusterinfo[(sv_type,scf,i)].append((p_val,nr_obs,mean_obs,window_size))
				break

		if new_cluster:
			self.clusters[(sv_type_observed, scf, self.index)] = {'max':pos,'min':pos,'pos':[pos]}
			self.clusterinfo[(sv_type_observed, scf, self.index)] = [(p_val,nr_obs,mean_obs,window_size)]
			self.index += 1


		if len(self.clusters) == 0:
			self.clusters[(sv_type_observed, scf, self.index)] = {'max':pos,'min':pos,'pos':[pos]}
			self.clusterinfo[(sv_type_observed, scf, self.index)] = [(p_val,nr_obs,mean_obs,window_size)]
			self.index += 1

	def get_final_bp_info(self):
		self.final_bps = []
		for region in self.clusters:

			avg_window_size = sum(map(lambda x:x[-1], self.clusterinfo[region]))/len(self.clusterinfo[region])
			start_pos = self.clusters[region]['min'] #- window_size
			end_pos = self.clusters[region]['max'] + self.clusterinfo[region][-1][-1]

			#n = len(self.clusters[region])
			#median_basepair = self.clusters[region][n/2]
			region_pvalues = map(lambda x:x[0], self.clusterinfo[region])
			avg_region_pvalue = sum(region_pvalues)/len(region_pvalues)
			region_nr_obs = map(lambda x:x[1], self.clusterinfo[region])
			avg_region_nr_obs = sum(region_nr_obs)/len(region_nr_obs)
			region_mean_obs = map(lambda x:x[2], self.clusterinfo[region])
			avg_region_mean_obs = sum(region_mean_obs)/len(region_mean_obs)

			self.final_bps.append( (region[1], 'GetDistr', 'FCD', start_pos, end_pos, avg_region_pvalue,'.','.','type:{0};avg_nr_span_obs:{1};mean_obs_isize:{2};window_size:{3}'.format(region[0], avg_region_nr_obs, avg_region_mean_obs, avg_window_size) ) )
			

	def __str__(self):
		"""
			Prints GFF/GTF file of expansion and contraction regions breakpoint regions
		"""
		output_string= 'seqname\tsource\tfeature\tstart\tend\tscore(p_val)\tstrand\tframe\tattribute\n'
		
		for seqname, source, feature, start, end, avg_p_val, strand, frame, attribute in self.final_bps:
			if start > self.param.mean and end < self.param.scaf_lengths[seqname] - self.param.mean:
				output_string += '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format(seqname, source, feature, int(start), int(end), avg_p_val, strand, frame, attribute) 
		return output_string

def AdjustInsertsizeDist(mean_insert, std_dev_insert, insert_list):
    k = MaxObsDistr(len(insert_list), 0.95)
    filtered_list = list(filter((lambda x : (x < mean_insert + k * std_dev_insert and x > mean_insert - k * std_dev_insert)), insert_list))
    if len(insert_list) > len(filtered_list):
        return(True, filtered_list)
    else:
        return(False, filtered_list)


def calc_p_values(bamfile,outfile,param, info_file,assembly_dict):

	p_values = []
	with pysam.Samfile(bamfile, 'rb') as bam:

		#sample true distribution
		param.sample_distribution(bam, info_file)

		# start scoring
		#reference_tids = map(lambda x: bam.gettid(x),bam.references )
		reference_lengths = map(lambda x: int(x), bam.lengths)
		scaf_dict = dict(zip(bam.references, reference_lengths))
		bam_filtered = ifilter(lambda r: r.flag <= 200, bam)
		current_scaf = -1
		print >> info_file, scaf_dict
		param.scaf_lengths = scaf_dict
		prev_coord = (0,0)
		duplicate_count = 0
		for i,read in enumerate(bam_filtered):
			if read.is_unmapped:
				continue
			else:
				current_coord = read.pos
				current_ref = bam.getrname(read.tid)

			coord1 = read.pos
			coord2 = read.mpos
			if (coord1, coord2) == prev_coord:
				continue
			else:
				prev_coord = (coord1, coord2)


			if (i + 1) %100000 == 0:
				# print i
				print >> info_file, '#Processing read:{0}'.format(current_coord)
				
			if (i+1) % 100000 == 0:
				print >> info_file, '#removed {0} duplicates, processed {1} reads.'.format(duplicate_count,i)
			# # # 	print 'extra!'
			# 	break

			# initialize read container for new scaffold
			if current_ref != current_scaf:
				print >> info_file, current_ref
				container = []
				scaf_length = scaf_dict[current_ref]

				# store positions in reverse to reduce complexity when removing items from lits in
				# python. We remove the last item and so on
				for i in range(scaf_length,0,-1):
					container.append(ReadContainer(i))
				current_scaf = current_ref 
			# the read pairs we want to use for calculating FCD
			# also, we follow reaprs suggestion and remove read pairs that are further away than
			# 1.5*param.mean from the position of interest. First we can safely remove anything
			# bigger than 3*param.mean (because at least one read of this read pair
			# is going to be further away than 1.5*mean from the position of interest in this read pair )
			if is_proper_aligned_unique_innie(read) and abs(read.tlen) <= 3*param.mean:
				#if abs(read.tlen) > 700:
				#	print read.tlen, read.pos,current_scaf
				#if current_scaf == 'scf_gap0_errorsize75' and (read.pos > 2352 or read.mpos > 2350):
				#	print current_scaf, read.pos, read.mpos
				if read.aend >= scaf_length or read.aend < 0 or read.mpos +read.rlen > scaf_length or read.pos < 0:
					#print 'Read coordinates outside scaffold length for {0}:'.format(current_scaf), read.aend, read.aend, read.mpos +read.rlen, read.pos 
					continue
	
				# Here we only add the observations that are not
				# further away than the pairs 1.5*param.mean < obs < 3*param.mean 

				# if abs(read.tlen) > 1.5*param.mean:
				# 	excluded_region_size = int(abs(read.tlen) - 1.5*param.mean)
				# 	#print 'LOOOOOOL'
				# else:
				# 	excluded_region_size = 0

				excluded_region_size = 0
				
				if read.tlen > 0:
					inner_start_pos = read.aend
					inner_end_pos = read.mpos
					for pos in range(inner_start_pos + excluded_region_size, inner_end_pos - excluded_region_size):
						try:
							container[scaf_length - pos].add_read(read)
						#container[pos].add_read(read2)
						except IndexError:
							pass
							#print 'Tried adding read pair to position {0}. On scaffold {1} with length {2}, with container of size {3}'.format(pos, current_scaf, scaf_length, len(container)) 
				else:

					inner_start_pos = read.mpos +read.rlen
					inner_end_pos = read.pos
					for pos in range(inner_start_pos + excluded_region_size, inner_end_pos - excluded_region_size):
						try:
							container[scaf_length - pos].add_read(read)
						except IndexError:
							pass
							#print 'Tried adding read pair to position {0}. On scaffold {1} with length {2}, with container of size {3}'.format(pos, current_scaf, scaf_length, len(container)) 

						#container[pos].add_read(read2)

			# write positions out to file
			if  current_coord > scaf_length - len(container):

				while scaf_length - len(container) < current_coord:

					# get true distribution
					if container[-1].position % 1000 == 0:
						print 'position', container[-1].position
					sequence_in_window = assembly_dict[ current_ref ][container[-1].position - int(1.5*param.mean) : container[-1].position + int(1.5*param.mean) ]
					#p = re.compile("[Nn]+")
					gap_coordinates = []
					#for m in p.finditer(sequence_in_window):
					#	gap_coordinates.append((m.start() - int(1.5*param.mean) ,m.end() - int(1.5*param.mean) ))

					true_distribution = param.get_correct_ECDF(outfile, gap_coordinates)
					container[-1].calc_observed_insert()
					KS_statistic, two_side_p_val = container[-1].calc_ks_test(true_distribution) 


					# do ks_2_sample
					if two_side_p_val > 0:
						p_values.append(two_side_p_val)
					# write chromosome, coord, p_val to file
					container[-1].write_pval_to_file(outfile,current_ref)
					del container[-1]



class Window(object):
	"""docstring for Window"""
	def __init__(self,pval_thresh, max_size):
		super(Window, self).__init__()
		self.queue = deque() 
		self.avg_mean = None
		self.nr_significant = 0
		self.nr_in_window = 0
		self.pval_thresh = pval_thresh
		self.avg_pval = 0
		self.max_window_size = max_size

	def update(self, pval,mean_isize):
		# if 0 <self.nr_significant < 100:
		# 	print self.nr_significant 
		# check = 0
		# if pval < 1.73242548129e-29:
		# 	check = 1
		# 	print 'HERE',self.avg_mean,self.nr_significant,self.nr_in_window,self.pval_thresh, self.avg_pval,self.max_window_size
		
		# initialize window

		if self.avg_mean == None and mean_isize > 0:
			# if check:
			# 	print 'here?', self.avg_mean,self.nr_significant,self.nr_in_window,self.pval_thresh, self.avg_pval,self.max_window_size
			self.avg_mean = mean_isize
			self.avg_pval = pval
			self.nr_in_window = 1
			if 0 <= pval < self.pval_thresh: 
				self.nr_significant = 1 
			else:
				self.nr_significant = 0
			self.queue.append((pval,mean_isize))
			#print 'looo'

		# update with new value

		if 0 <= pval < self.pval_thresh: 
			self.nr_significant += 1 
		self.avg_mean = (self.nr_in_window * self.avg_mean + mean_isize)/ float(self.nr_in_window+1)
		self.avg_pval = (self.nr_in_window * self.avg_pval + pval)/ float(self.nr_in_window+1)
		self.queue.append((pval,mean_isize))
		self.nr_in_window += 1

		# window is full

		if self.nr_in_window >= self.avg_mean or  self.nr_in_window >= self.max_window_size:
			# if check:
			# 	print 'be here!', self.avg_mean,self.nr_significant,self.nr_in_window,self.pval_thresh, self.avg_pval,self.max_window_size

			if self.is_significant():
				# if check:
				# 	print 'please here!!', self.avg_mean,self.nr_significant,self.nr_in_window,self.pval_thresh, self.avg_pval,self.max_window_size

				pval_left, mean = self.queue.popleft()
				# if check:
				# 	print pval_left

				self.avg_mean = (self.nr_in_window * self.avg_mean - mean)/ float(self.nr_in_window-1)
				self.avg_pval = (self.nr_in_window * self.avg_pval - pval_left)/ float(self.nr_in_window-1)
				self.nr_in_window -= 1
				if 0 <= pval_left < self.pval_thresh:
					self.nr_significant -=1
		
				return True
			else:
				pval_left, mean = self.queue.popleft()
				# if check:
				# 	print pval_left

				self.avg_mean = (self.nr_in_window * self.avg_mean - mean)/ float(self.nr_in_window-1)
				self.avg_pval = (self.nr_in_window * self.avg_pval - pval_left)/ float(self.nr_in_window-1)
				self.nr_in_window -= 1
				if 0 <= pval_left < self.pval_thresh:
					self.nr_significant -=1		
				return False		



	def is_significant(self):
		"""
		The window has a misassembly according to reapers thresholds, that is 
		, at least 80% of the bases in the window has a p-value under pval. Outfolder
		window size is defined as min(reaprs' window size, avg fragment length in window)
		This is because a contraction could not be found otherwise since it will reduce 
		the all insert sizes spanning over.
		"""
		# window is large, need to drop leftmost position and 
		# check if window is significant
	
		if self.nr_significant /float(self.nr_in_window) >= 0.8:
			return True
		else:
			return False



		
def get_misassemly_regions(pval_file, param, info_file):
	#information = [ map(lambda x: float(x), line.strip().split()[1:] ) for line in pval_file.readlines()[:10]]
	current_seq = -1
	sv_container = BreakPointContainer(param)
	for line in pval_file.readlines():
		scf, pos, pos_p_val, n_obs, mean, stddev = line.strip().split()
		if float(pos_p_val) == -1 or float(mean) < 0:
			current_seq = -1
			continue

		if (scf != current_seq and pos >= param.max_window_size):
			current_seq = scf
			window = Window(param.pval, param.max_window_size)
			window.update(float(pos_p_val), float(mean))
		
		else:
			became_significant = window.update(float(pos_p_val), float(mean))
			if became_significant:
				#print 'Sign'
				if window.avg_mean > param.adjusted_mean:
					sv_container.add_bp_to_cluster(current_seq, int(pos), window.avg_pval, int(n_obs), window.avg_mean, 'expansion', min(window.nr_in_window, window.max_window_size))
				else:
					sv_container.add_bp_to_cluster(current_seq, int(pos), window.avg_pval, int(n_obs), window.avg_mean, 'contraction', min(window.nr_in_window, window.max_window_size))

	print 'read pval file'

	return sv_container

def scan_bam(bam_file, assembly_file, outfolder):
	param = Parameters()
	if not os.path.exists(args.outfolder):
		os.makedirs(args.outfolder)
	info_file = open(os.path.join(args.outfolder,'info.txt'),'w')
	pval_file_out = open(os.path.join(args.outfolder,'p_values.txt'),'w')
	assembly_dict = ReadInContigseqs(open(assembly_file,'r'),param.window_size)

	calc_p_values(bam_file, pval_file_out, param, info_file,assembly_dict)
	pval_file_out.close()


def cluster_pvals(pval_file_in, outfile, info_file ,assembly_file,p_val_threshold, window_size):
	param = Parameters()
	vals = filter( lambda line: line[0] != '#', open(info_file,'r').readlines())[0:2]
	print vals
	[mean,stddev] =  vals[0].strip().split()
	[adjusted_mean, adjusted_stddev] = vals[1].strip().split()
	param.mean, param.stddev, param.adjusted_mean,param.adjusted_stddev = float(mean), float(stddev), float(adjusted_mean), float(adjusted_stddev)
	param.pval = p_val_threshold
	print param.mean, param.stddev, param.adjusted_mean,param.adjusted_stddev
	# if window_size >= 1000:
	# 	param.max_window_size = window_size/2 
	# else:
	# 	param.max_window_size = window_size
	param.max_window_size = window_size
	gff_file = open(outfile,'w')
	pval_file_in = open(pval_file_in,'r')
	assembly_dict = ReadInContigseqs(open(assembly_file,'r'), param.max_window_size)

	param.scaf_lengths = dict(map( lambda (key,val): (key,len(val)),assembly_dict.iteritems()))
	print param.scaf_lengths
	sv_container =  get_misassemly_regions(pval_file_in,param, info_file) #open(args.pval_file,'r')
	sv_container.get_final_bp_info()
	print >> gff_file, str(sv_container)


def main_pipline(args):
		scan_bam(args.bampath, args.assembly_file, args.outfolder)
		cluster(args.pval_file, args.outfile, args.info_file ,args.pval, args.window_size)





	# param.pval = args.pval


	
	# param.mean = 2570.0
	# param.adjusted_mean = 3278.6




if __name__ == '__main__':

	# create the top-level parser
	parser = argparse.ArgumentParser()#prog="Infer variants with simple p-value test using theory of GetDistr - proof of concept.")
	#parser.add_argument('--foo', action='store_true', help='help for foo arg.')
	subparsers = parser.add_subparsers(help='help for subcommand')

	# create the parser for the "pipeline" command
	pipeline = subparsers.add_parser('pipeline', help='Run the entire pipeline')
	pipeline.add_argument('a', type=str, help='help for bar, pipeline')

	pipeline.add_argument('bampath', type=str, help='bam file with mapped reads. ')
	pipeline.add_argument('assembly_file', type=str, help='Fasta file with assembly/genome. ')
	pipeline.add_argument('outfolder', type=str, help='Outfolder. ')
	pipeline.add_argument('window_size', type=int, help='Window size ')
	pipeline.add_argument('outfolder', type=str, help='Outfolder. ')
	pipeline.add_argument('pval', type=float, help='p-value threshold for calling a variant. ')
	pipeline.set_defaults(which='pipeline')
	
	# create the parser for the "scan bam" command
	scan_bam = subparsers.add_parser('scan_bam', help='Scan bam file and calculate pvalues for each base pair')
	scan_bam.add_argument('bampath', type=str, help='bam file with mapped reads. ')
	scan_bam.add_argument('assembly_file', type=str, help='Fasta file with assembly/genome. ')
	scan_bam.add_argument('outfolder', type=str, help='Outfolder. ')
	scan_bam.set_defaults(which='scan_bam')


	# create the parser for the "cluster" command
	cluster = subparsers.add_parser('cluster_pvals', help='Takes a pvalue file and clusters them into significan regions')
	cluster.add_argument('pval_file', type=str, help='p-val file from "scan_bam" output. ')
	cluster.add_argument('info_file', type=str, help='info file from "scan_bam" output. ')
	cluster.add_argument('assembly_file', type=str, help='Fasta file with assembly/genome. ')
	cluster.add_argument('pval', type=float, help='p-value threshold for calling a variant. ')
	cluster.add_argument('window_size', type=int, help='Window size ')
	cluster.add_argument('outfile', type=str, help='Outfile. ')
	cluster.set_defaults(which='cluster')

	# parser.add_argument('mean', type=int, help='mean insert size. ')
	# parser.add_argument('stddev', type=int, help='Standard deviation of insert size ')

	
	args = parser.parse_args()
	#print args

	if args.which == 'pipeline':
		main_pipline(args)
	elif args.which == 'scan_bam':
		scan_bam(scan_bam.bampath, scan_bam.assembly_file, scan_bam.outfolder)
	elif args.which == 'cluster':
		cluster_pvals(args.pval_file, args.outfile, args.info_file, args.assembly_file ,args.pval, args.window_size)


	#main(args)



        
