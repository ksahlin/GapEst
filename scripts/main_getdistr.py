import sys,os,subprocess

import argparse

def main(args):
    plotfolder = '/tmp/gapest_simulated/plots'
    if not os.path.exists(plotfolder):
        os.makedirs(plotfolder)

    for distr in ['normal', 'uniform', 'mix']:
        ## simulate instance
        #os.popen("python scripts/simulate_data.py 2000000 1000 4000 {0} 50 100 /tmp/gapest_simulated/{0}/ -sort  --mean 2000 --sd 500 --contigs /tmp/gapest_simulated/ctgs.fa --genome /tmp/gapest_simulated/genome.fa ".format(distr))
        ## Get true gaps
        #os.popen("python scripts/evaluate_gapest.py --getgaps /tmp/gapest_simulated/genome.fa /tmp/gapest_simulated/ctgs.fa /tmp/gapest_simulated/{0}/true_gaps".format(distr))
        ## Get estimated gaps     
        #for est_type in ['gapest', 'naive']:
        for est_type in [ 'naive']:
            if not os.path.exists("/tmp/gapest_simulated/{0}/{1}/".format(distr,est_type)):
                os.makedirs("/tmp/gapest_simulated/{0}/{1}/".format(distr,est_type))
            if est_type == 'naive':
                os.popen("python /Users/ksahlin/Documents/workspace/GapEst/src/Main.py 1 -c /tmp/gapest_simulated/ctgs.fa -f /tmp/gapest_simulated/{0}/mapped.bam -m 2000 -s 500 -e 10 -r 30 --naive 1 > /tmp/gapest_simulated/{0}/{1}/{1}.gaps ".format(distr,est_type) )
            else:
                os.popen("python /Users/ksahlin/Documents/workspace/GapEst/src/Main.py 1 -c /tmp/gapest_simulated/ctgs.fa -f /tmp/gapest_simulated/{0}/mapped.bam -m 2000 -s 500 -e 10 -r 30 >  /tmp/gapest_simulated/{0}/{1}/{1}.gaps".format(distr,est_type) )

            ## plot results
            os.popen("python scripts/evaluate_gapest.py --comparegaps /tmp/gapest_simulated/{0}/true_gaps/truegaps.gaps /tmp/gapest_simulated/{0}/{1}/{1}.gaps {0}_{1} {2}".format(distr,est_type, plotfolder))
	

if __name__ == '__main__':
    ##
    # Take care of input
    parser = argparse.ArgumentParser(description="Simulate paired end reads with a given mean and standard deviation on mean insert size. A given coverage should also be provided.")
    #parser.add_argument('distr', type=str, help='Distribution of insert sizes. ')
    #parser.add_argument('coverage', type=int, help='Coverage. ')
    #parser.add_argument('read_length', type=int, help='Length of read fragments within a pair. ')
    #parser.add_argument('output_path', type=str, help='path to folder output. ')
    #parser.add_argument( '-sort', dest='sort', action='store_true', default=False, help='Coordinate sort the reads in the bam file' )
    #parser.add_argument( '-sam', dest='sam', action='store_true', default=False, help='Output a samfile (default is bam)' )
    #parser.add_argument('--threads', type=str, dest='threads', default='8', required=False, help='Number of threads for bwa mem.')

    #parser.add_argument('--contigs', dest='contigs', type=str, default=False, help='Path to already generated contigs. ')
    #parser.add_argument('--genome', dest='genome', type=str, default=False, help='Path to already generated genome. ')

    #parser.add_argument('--mean', dest='mean',type=int, default=None, help='Mean insert. ')
    #parser.add_argument('--sd', dest='sd', type=int, default=None, help='Stddev insert.')
    #parser.add_argument('--min_size', dest='min_size', type=int, default=None, help='Min insert size (only effective if distr=uniform). ')
    #parser.add_argument('--max_size', dest='max_size', type=int, default=None, help='Max insert size (only effective if distr=uniform). ')
   

    args = parser.parse_args()


    main(args)








