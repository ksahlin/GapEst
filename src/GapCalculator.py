import sys
import random 

import numpy as np
import Contig,Scaffold
from scipy.special import erf
from scipy.stats import norm
from scipy.constants import pi
from math import exp

from mathstats.normaldist.normal import MaxObsDistr
from mathstats.normaldist.truncatedskewed import param_est

def subsample_gap_est(observations,mean,sigma, effective_read_length, c1_len, c2_len):
    
    gap_estimations = []
    #print observations, c1_len,c2_len
    for i in range(20):
        sample_size = max(int(len(observations)*random.uniform(0,1)),1) # let sample size be 20% or the initial nr of observations
        subsample = [ observations[i] for i in sorted(random.sample(xrange(len(observations)),sample_size )) ]
        o_mean = int(sum(subsample)/len(subsample))
        d_hat = int(param_est.GapEstimator(mean,sigma, effective_read_length, o_mean, c1_len, c2_len))
        gap_estimations.append(d_hat)

    gap_naive = mean - sum(observations)/len(observations)
    index = np.argmin(map(lambda x: abs(gap_naive-x), gap_estimations))
    # avg_gap = sum(gap_estimations) / len(gap_estimations)
    #print gap_estimations
    return(gap_estimations[index])


def AdjustInsertsizeDist(mean_insert, std_dev_insert, insert_list):
    k = MaxObsDistr(len(insert_list), 0.95)
    filtered_list = list(filter((lambda x : (x < mean_insert + k * std_dev_insert and x > mean_insert - k * std_dev_insert)), insert_list))
    if len(insert_list) > len(filtered_list):
        return(True, filtered_list)
    else:
        return(False, filtered_list)

#def has_mismappings(sorted_observations):
#    max_dist_between_obs = max(map(lambda x,y: y-x, izip(sorted_observations[:-1], sorted_observations[1:]) ) )

def remove_misalignments(sorted_observations,edge_support):
    if len(sorted_observations) < edge_support:
            return sorted_observations
    n_isize = float(len(sorted_observations))
    mean_isize = sum(sorted_observations)/n_isize
    std_dev_isize =  (sum(list(map((lambda x: x ** 2 - 2 * x * mean_isize + mean_isize ** 2), sorted_observations))) / (n_isize - 1)) ** 0.5
    #print '#Mean before filtering :', mean_isize
    #print '#Stddev before filtering: ', std_dev_isize
    extreme_obs_occur = True
    while extreme_obs_occur:
        extreme_obs_occur, filtered_list = AdjustInsertsizeDist(mean_isize, std_dev_isize, sorted_observations)
        if len(filtered_list)< edge_support:
            return filtered_list
        n_isize = float(len(filtered_list))
        mean_isize = sum(filtered_list) / n_isize
        std_dev_isize = (sum(list(map((lambda x: x ** 2 - 2 * x * mean_isize + mean_isize ** 2), filtered_list))) / (n_isize - 1)) ** 0.5
        sorted_observations = filtered_list

    #print '#Mean converged:', mean_isize
    #print '#Std_est converged: ', std_dev_isize
    return sorted_observations

def GapEstimator(G,Contigs,Scaffolds,mean,sigma,effective_read_length,edge_support,bayesian, naive, reinforced):
    ## effective_read_length is here that average effective_read_length minus max allowed softclipps (i.e. the smallest part that can be aligned)
    gap_mean=0
    gap_obs=0
    gap_obs2=0
    gap_counter=0
    print 'CONTIG1\tCONTIG2\tGAP_ESTIMATION\tC1DIR\tC2DIR\tNUMBER_OF_OBSERVATIONS\tC1LEN\tC2LEN\tWARNINGS/ERRORS'
    for edge in G.edges_iter():
        if G[edge[0]][edge[1]]['nr_links'] != None:
            c1=edge[0][0]
            c2=edge[1][0]

            c1_dir = edge[0][1]
            c2_dir = edge[1][1]
            #print Scaffolds[c1].contigs[0].name,Scaffolds[c2].contigs[0].name

            if c1_dir == 'R':
                c1_orientation = '+'
            else:
                c1_orientation = '-'
            if c2_dir == 'R':
                c2_orientation = '+'
            else:
                c2_orientation = '-'

            c1_len=G.node[edge[0]]['length']
            c2_len=G.node[edge[1]]['length']
            obs_list=G.edge[edge[0]][edge[1]]['gap_dist']
            nr_links=G.edge[edge[0]][edge[1]]['nr_links']
            min_obs_pos = min(map( lambda x: x[0], G.edge[edge[0]][edge[1]]['obs_pos']))
            #pre check for large deviations in obs list
            sorted_observations = sorted(obs_list)
            filtered_observations = remove_misalignments(sorted_observations,edge_support)
            smallest_obs_mean = sum(sorted_observations[0:10])/10.0
            largest_obs_mean = sum(sorted_observations[-10:])/10.0
            #print largest_obs_mean,smallest_obs_mean

            if (len(filtered_observations) >= edge_support) and (largest_obs_mean-smallest_obs_mean < 6*sigma):
                if bayesian: 
                    d_hat,stdErr=CalcMLvaluesOfdGeneralBayesian(filtered_observations,mean,sigma,effective_read_length,c1_len,c2_len,nr_links)
                
                elif naive:
                    d_hat,stdErr = mean - int(sum(filtered_observations)/len(filtered_observations)), 0 
                elif reinforced:
                    n_isize = len(filtered_observations)
                    o_mean = int(sum(filtered_observations)/n_isize)
                    o_stddev = (sum(list(map((lambda x: x ** 2 - 2 * x * o_mean + o_mean ** 2), filtered_observations))) / (n_isize - 1)) ** 0.5

                    d_hat,stddev_est = param_est.GapEstimator_reinforced(mean,sigma, effective_read_length, o_mean, o_stddev, c1_len, c2_len)
                    d_hat = int(d_hat)
                    stddev_est = int(stddev_est)
                    # print 'ESTIMATED CONDITIONAL STDDEV: {0}, OBSERVED:{1}'.format( stddev_est, o_stddev)
                    # if stddev_est > o_stddev:
                    #     print 'Gap est:{0} -> stddev adjust would give higher gap estimation'.format(d_hat)
                    # else:
                    #     print 'Gap est:{0} -> stddev adjust would give lower gap estimation'.format(d_hat)
                else:
                    o_mean = int(sum(filtered_observations)/len(filtered_observations))
                    d_hat = int(param_est.GapEstimator(mean,sigma, effective_read_length, o_mean, c1_len, c2_len))
                    
                    if d_hat > mean+3*sigma -2*effective_read_length :
                        #get average gap from subsampling
                        d_hat = subsample_gap_est(filtered_observations, mean, sigma, effective_read_length, c1_len, c2_len)
                    # no negative observations, we cannot have two contigs 
                    # overlapping more that the smallest observation
                    if d_hat <0 and -d_hat > min_obs_pos: 
                        d_hat = 2*effective_read_length - min_obs_pos
                    # smallest_obs = min(filtered_observations)
                    # # no negative observations
                    # if smallest_obs +  d_hat < 2*effective_read_length:
                    #     d_hat = 2*effective_read_length - smallest_obs

                    #d_hat,stdErr=CalcMLvaluesOfdGeneral(filtered_observations,mean,sigma,effective_read_length,c1_len,c2_len,nr_links) 
              
                
                gap_obs+=d_hat
                gap_obs2+=d_hat**2
                gap_counter+=1

                if c1_len < 2*sigma and c2_len < 2*sigma:
                    print Scaffolds[c1].contigs[0].name + '\t'+ Scaffolds[c2].contigs[0].name+ '\t'+str(d_hat) + '\t'+str(c1_orientation)+'\t'+str(c2_orientation)+ '\t'+str(nr_links)+'\t'+str(c1_len)+'\t'+str(c2_len)+'\tw1,c1_len={0}c2_len={1},{2}'.format(c1_len,c2_len,filtered_observations)
                elif c1_len < (2*sigma + effective_read_length) or c2_len < (2*sigma + effective_read_length):
                    print Scaffolds[c1].contigs[0].name + '\t'+ Scaffolds[c2].contigs[0].name+ '\t'+str(d_hat) + '\t'+str(c1_orientation)+'\t'+str(c2_orientation)+ '\t'+str(len(filtered_observations))+'\t'+str(c1_len)+'\t'+str(c2_len)+'\tw3:c1_len={0}c2_len={1},{2}'.format(c1_len, c2_len, len(filtered_observations))    
                elif filtered_observations < sorted_observations:
                    print Scaffolds[c1].contigs[0].name + '\t'+ Scaffolds[c2].contigs[0].name+ '\t'+str(d_hat) + '\t'+str(c1_orientation)+'\t'+str(c2_orientation)+ '\t'+str(len(filtered_observations))+'\t'+str(c1_len)+'\t'+str(c2_len)+'\tw2:prev_nr_links:{0},links_after_filter:{1},c1_len={3}c2_len={4},{2}'.format(nr_links,len(filtered_observations), filtered_observations,c1_len,c2_len)                
                else:
                    print Scaffolds[c1].contigs[0].name + '\t'+ Scaffolds[c2].contigs[0].name+ '\t'+str(d_hat) + '\t'+str(c1_orientation)+'\t'+str(c2_orientation)+ '\t'+str(nr_links)+'\t'+str(c1_len)+'\t'+str(c2_len)+'\t-,c1_len={0}c2_len={1},{2}'.format(c1_len,c2_len,filtered_observations)


            elif nr_links < edge_support:
                print Scaffolds[c1].contigs[0].name + '\t'+ Scaffolds[c2].contigs[0].name+ '\t.' + '\t'+str(c1_orientation)+'\t'+str(c2_orientation)+ '\t.\te1'
            elif largest_obs_mean-smallest_obs_mean < 6*sigma:
                print Scaffolds[c1].contigs[0].name + '\t'+ Scaffolds[c2].contigs[0].name+ '\t.' + '\t'+str(c1_orientation)+'\t'+str(c2_orientation)+ '\t.\te2'
                    
    print 'w1 : Both contig lengths were smaller than 2*std_dev of lib (heuristic threshold set by me from experience). This can give shaky estimations in some cases.'
    print 'w2 : GapEst filtered out at least one extreme outlier in observation (potential mismapping or misassembly) before calculating gap esitmation.'
    print 'w3 : One contig is smaller than sigma + effective_read_length. Estimation can be shaky. Use naive estimation instead'
    print 'e1 : No gap was calculated. Number of links were lower than specified min nr of links parameter: -e <min nr links> (default 10). Lower this value if necessary (estimations may become unstable)'
    print 'e2 : No gap was calculated. The spread of the links throughout the contig is to far (i.e largest_obs_mean-smallest_obs_mean < 6*sigma ), suggesting errors in mapping on this region.',
                
    gap_mean=gap_obs/float(gap_counter)
    gap_sigma=gap_obs2-gap_counter*gap_mean**2

def funcDGeneral(obs_list,d,mean,stdDev,c1Len,c2Len,readLen):
    #get observation    
    #data_observation=(nr_links*mean -int(sum(obs_list)))/float(nr_links)
    c_min=min(c1Len,c2Len)
    c_max=max(c1Len,c2Len)

    def Nominatortest(d,c_min,c_max,c1Len,c2Len,readLen):
        nomin1=( erf((c_min+c_max+d+1-mean)/(2**0.5*float(stdDev))) +erf((mean-d-c_max-readLen-1)/(2**0.5*float(stdDev))))*(pi/2)**0.5
        nomin2=-(erf((c_min+d+readLen+1-mean)/(2**0.5*float(stdDev)))+erf((mean-d-2*readLen+1)/(2**0.5*float(stdDev))))*(pi/2)**0.5
        nomin=nomin1+nomin2
        case_nomin=0
        return nomin,case_nomin
    
    def Denominator(d,c1Len,c2Len,readLen):
        #term 1,2 and 3 denodes what part of the function we are integrating term1 for first (ascending), etc...
        term2=(c_min-readLen+1)/2.0*(erf((c_max+d+readLen-mean)/((2**0.5)*stdDev))- erf((c_min+d+readLen-mean)/((2**0.5)*stdDev))   )
    
        first=-((pi/2)**0.5)*(d+2*readLen-mean-1)*( erf((c_min+d+readLen-mean)/(2**0.5*float(stdDev))) - erf((d+2*readLen-1-mean)/(2**0.5*float(stdDev)))  )
        second=stdDev*( exp(-( (d+2*readLen-1-mean)**2)/(float(2*stdDev**2))) - exp(-( (c_min+d+readLen-mean)**2)/(float(2*stdDev**2)))) 
        term1=first+second

        first=((pi/2)**0.5)*(c_min+c_max+d-mean+1)*( erf((c_min+c_max+d-mean)/(2**0.5*float(stdDev))) - erf((c_max+readLen+d-mean)/(2**0.5*float(stdDev)))  )
        #print 'First: ',first
        second=stdDev*( exp(-( (c_min+c_max+d-mean)**2)/(float(2*stdDev**2))) - exp(-( (c_max+readLen+d-mean)**2)/(float(2*stdDev**2))))
        #print 'Second: ',second
        term3=first+second
        denom=term1+term2+term3
        #print term1,term2,term3
        return denom
    
    denominator=Denominator(d,c1Len,c2Len,readLen)
    nominator,case_nomin=Nominatortest(d,c_min,c_max,c1Len,c2Len,readLen)
    Aofd=nominator/denominator
    func_of_d=d+Aofd*stdDev**2
    return func_of_d

def CalcMLvaluesOfdGeneral(obs_list,mean,stdDev,readLen,c1Len,c2Len,nr_links):
    #get observation    
    data_observation=(nr_links*mean -int(sum(obs_list)))/float(nr_links)
    #do binary search among values
    d_upper= max(int(mean+4*stdDev-2*readLen), 2*readLen)
    d_lower= min(-int(mean+4*stdDev-2*readLen), -2*readLen) 
    while d_upper-d_lower>1:
        d_ML=(d_upper+d_lower)/2.0
        func_of_d=funcDGeneral(obs_list,d_ML,mean,stdDev,c1Len,c2Len,readLen)
        #print d_ML,func_of_d
        if func_of_d>data_observation:
            d_upper=d_ML
        else:
            d_lower=d_ML

    d_ML=(d_upper+d_lower)/2.0
    return int(round(d_ML,0)),0

def CalcMLvaluesOfdGeneralBayesian(obs_list,mean,stdDev,readLen,c1Len,c2Len,nr_links):
    #get observation    
    data_observation=(nr_links*mean -int(sum(obs_list)))/float(nr_links)
    #do binary search among values
    d_upper=int(mean+10*stdDev-2*readLen)
    d_lower=-10*stdDev
    while d_upper-d_lower>1:
        d_MAP=(d_upper+d_lower)/2.0
        func_of_d=BayesianEquation(obs_list,d_MAP,mean,stdDev,c1Len,c2Len,readLen)
        #print d_MAP,func_of_d
        if func_of_d>data_observation:
            d_upper=d_MAP
        else:
            d_lower=d_MAP

    d_MAP=(d_upper+d_lower)/2.0
    return int(round(d_MAP,0)),0  

def BayesianEquation(obs_list,d,mean,stdDev,c1Len,c2Len,readLen):
    ML_value = funcDGeneral(obs_list,d,mean,stdDev,c1Len,c2Len,readLen)
    prior_value = stdDev**2 * (1/ (1- norm.cdf(d,mean,stdDev))) * norm.pdf(d,mean,stdDev)

    return ML_value + prior_value
 
