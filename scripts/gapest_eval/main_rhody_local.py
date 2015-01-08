import sys,os,subprocess

import argparse

def main(args):
    plotfolder = '/tmp/rhody_gapest_eval/plots'
    main_out = '/tmp/rhody_gapest_eval/'
    bam_folder = '/Users/ksahlin/Downloads/aligned/rhody/'
    ctg_folder = '/Users/ksahlin/Work/KTH/BESST-scaffolder/evaluation/Rodhobacter/Assembly/'
    genome = '/Users/ksahlin/Work/KTH/BESST-scaffolder/evaluation/Rodhobacter/reference.fasta'
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
        #os.popen("python /Users/ksahlin/Documents/workspace/GapEst/scripts/evaluate_gapest.py --getgaps {0} {1} {2}/true_gaps 7000".format(genome, ctgs, outdir))
        ## Get estimated gaps     
        for est_type in ['gapest', 'naive']:
            if est_type == 'naive':
                os.popen("python /Users/ksahlin/Documents/workspace/GapEst/src/Main.py 1 -c {0} -f {1} -m 2563 -s 1336 -e 5 -r 101 --naive 1 > {2} ".format(ctgs, bam_file, os.path.join(outdir, est_type+'.gaps') ))
            else:
                os.popen("python /Users/ksahlin/Documents/workspace/GapEst/src/Main.py 1 -c {0} -f {1} -m 2563 -s 1336 -e 5 -r 101 >  {2} ".format(ctgs,bam_file, os.path.join(outdir, est_type+'.gaps')) )

            ## plot results
            os.popen("python /Users/ksahlin/Documents/workspace/GapEst/scripts/evaluate_gapest.py --comparegaps {4}/true_gaps/truegaps.gaps {0} {1}_{2} {3}".format(os.path.join(outdir, est_type+'.gaps'), assembly, est_type, plotfolder,outdir))
	

if __name__ == '__main__':
    ##
    # Take care of input
    parser = argparse.ArgumentParser(description="Gapestimations of rhodobacter.")

    args = parser.parse_args()


    main(args)








