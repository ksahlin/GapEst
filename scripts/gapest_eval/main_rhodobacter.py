import sys,os,subprocess

import argparse

def main(args):
    plotfolder = '/home/kris/Work/GetDistr/gapest_eval/plots'
    main_out = '/home/kris/Work/GetDistr/gapest_eval/'
    bam_folder = '/proj/b2013169/private/data/bacgenomes/aligned/rhody/'
    ctg_folder = '/proj/b2013169/private/data/gage.cbcb.umd.edu/data/Rhodobacter_sphaeroides/Assembly/'
    genome = '/proj/b2013169/private/data/bacgenomes/references/modified_references/rhody.reference.fasta'
    if not os.path.exists(plotfolder):
        os.makedirs(plotfolder)

    for assembly in [ 'ABySS', 'ABySS2', 'Allpaths-LG', 'Bambus2', 'CABOG', 'MSR-CA', 'SGA', 'SOAPdenovo', 'Velvet']:
        ## Get true gaps
        tmp_path = os.path.join(ctg_folder,assembly)
        ctgs = os.path.join(tmp_path,'genome.ctg.fasta')
        bam_file = os.path.join(bam_folder, assembly+'.bam')
        outdir = os.path.join(main_out, assembly)

        if not os.path.exists(outdir):
            os.makedirs(outdir)        
        os.popen("python /home/kris/git_repos/GapEst/scripts/evaluate_gapest.py --getgaps {0} {1} {2}/true_gaps".format(genome, ctgs, outdir))
        ## Get estimated gaps     
        for est_type in ['gapest', 'naive']:
            #est_type_out = os.path.join(outdir, est_type)
            #if not os.path.exists(est_type_out):
            #    os.makedirs(est_type_out)
            if est_type == 'naive':
                os.popen("python /home/kris/git_repos/GapEst/src/Main.py 1 -c {0} -f {1} -m 2600 -s 1325 -e 5 -r 30 --naive 1 > {2} ".format(ctgs, bam_file, os.path.join(outdir, est_type+'.gaps') ))
            else:
                os.popen("python /home/kris/git_repos/GapEst/src/Main.py 1 -c {0} -f {1} -m 2600 -s 1325 -e 5 -r 30 >  {2} ".format(ctgs,bam_file, os.path.join(outdir, est_type+'.gaps')) )

            ## plot results
            os.popen("python /home/kris/git_repos/GapEst/scripts/evaluate_gapest.py --comparegaps {4}/true_gaps/truegaps.gaps {0} {1}_{2} {3}".format(os.path.join(outdir, est_type+'.gaps'), assembly, est_type, plotfolder,outdir))
	

if __name__ == '__main__':
    ##
    # Take care of input
    parser = argparse.ArgumentParser(description="Gapestimations of rhodobacter.")
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








