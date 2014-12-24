import sys,os,subprocess

import argparse
from genomics_tools.file_formats import bam,fasta
from genomics_tools.simulate import genome,contigs,reads
from mapping import align
from genomics_tools.file_formats.various_annotations import to_AGP,to_GFF

def simulate_instance(args):
    print 'Started simulating'
    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path)
    
    if not args.contigs:
        genome_path = os.path.join(args.output_path, 'genome.fa')
        contig_path = os.path.join(args.output_path, 'ctgs.fa')
    else:
        genome_path = args.genome
        contig_path = args.contigs

    read1_path = os.path.join(args.output_path, 'reads1.fa')
    read2_path = os.path.join(args.output_path, 'reads2.fa')
    bam_path = os.path.join(args.output_path, 'mapped')

    if not args.contigs:
        #genome
        print args.genomelen
        g = genome.Genome([0.25]*4,args.genomelen,'genome1')
        g.genome()
        print >> open(genome_path,'w'), g.genome_fasta_format()

        #contigs
        ctgs = open(contig_path,'w')
        for ctg in contigs.generate_contigs(g.sequence,args.min_contig, args.max_contig, 0,3000):
            ctgs.write(ctg)
    else:
        g = genome.Genome([0.25]*4,args.genomelen,'genome1')
        #print genome_path, args.genomelen
        for acc,seq in fasta.fasta_iter(open(genome_path,'r')):
            #print acc, seq
            g.sequence = seq
            g.accession = acc
    #ctgs.write('>ctg0\n{0}\n'.format(g.sequence[0:args.burnin]))
    #for i,x in enumerate(range(args.burnin,args.genomelen,(args.contiglen + args.gaplen))):
    #	ctgs.write('>ctg{0}\n{1}\n'.format(i+1,g.sequence[x:x+args.contiglen]))

    #reads
    if args.distr == 'normal':
        lib = reads.DNAseq(args.read_length ,args.coverage, distribution=args.distr, mean=args.mean,stddev=args.sd)
        lib.simulate_pe_reads(g)
    elif args.distr == 'uniform':
        lib = reads.DNAseq(args.read_length ,args.coverage, distribution=args.distr, min_size=args.min_size,max_size=args.max_size)
        lib.simulate_pe_reads(g)
    elif args.distr == 'mix':
        lib_part1 = reads.DNAseq(args.read_length ,args.coverage/2, distribution='normal', mean=args.mean,stddev=args.sd)
        lib_part1.simulate_pe_reads(g)
        lib_part2 = reads.DNAseq(args.read_length ,args.coverage/2, distribution='uniform', min_size=(args.mean - 3*args.sd),max_size=(args.mean + 3*args.sd))
        lib_part2.simulate_pe_reads(g)
        # concatenate the reads from each distribution
        lib = reads.DNAseq(args.read_length ,args.coverage, distribution=args.distr, mean=args.mean,stddev=args.sd)
        lib.reads = lib_part1.reads + lib_part2.reads


    reads1 = open(read1_path,'w')
    reads2 = open(read2_path,'w')
    i=0
    for read in lib.fasta_format():
        if i%2==0:
            reads1.write(read)
        else:
            reads2.write(read)
        i+=1

    print 'Started mapping'
    #mapping
    align.map_paired_reads(read1_path, read2_path, contig_path, bam_path, args)
    #align.bwa_mem(read1_path, read2_path, contig_path, bam_path, args)

def main(args):
    successful_experiments = 0
    while successful_experiments < 1: 
        try:
            simulate_instance(args)
        except subprocess.CalledProcessError:
            continue

        successful_experiments += 1
	

if __name__ == '__main__':
    ##
    # Take care of input
    parser = argparse.ArgumentParser(description="Simulate paired end reads with a given mean and standard deviation on mean insert size. A given coverage should also be provided.")
    parser.add_argument('genomelen', type=int, help='Length of contigs. ')    
    parser.add_argument('min_contig', type=int, help='Length of contigs. ')
    parser.add_argument('max_contig', type=int, help='Length of contigs. ')
    parser.add_argument('distr', type=str, help='Distribution of insert sizes. ')
    parser.add_argument('coverage', type=int, help='Coverage. ')
    parser.add_argument('read_length', type=int, help='Length of read fragments within a pair. ')
    parser.add_argument('output_path', type=str, help='path to folder output. ')
    parser.add_argument( '-sort', dest='sort', action='store_true', default=False, help='Coordinate sort the reads in the bam file' )
    parser.add_argument( '-sam', dest='sam', action='store_true', default=False, help='Output a samfile (default is bam)' )
    #parser.add_argument( '-scafs', dest='scaffolds', action='store_true', default=False, help='scaffolds are simulated instead of contigs' )
    #parser.add_argument( '-errors', dest='errorsize', type=int, nargs='+', default=False, help='gap distance error' )
    #parser.add_argument( '-burnin', dest='burnin', type=int, default=False, help='Estimation window' )
    #parser.add_argument( '-nrgaps', dest='nrgaps', type=int, default=False, help='Number of gaps' )
    parser.add_argument('--threads', type=str, dest='threads', default='8', required=False, help='Number of threads for bwa mem.')

    parser.add_argument('--contigs', dest='contigs', type=str, default=False, help='Path to already generated contigs. ')
    parser.add_argument('--genome', dest='genome', type=str, default=False, help='Path to already generated genome. ')

    parser.add_argument('--mean', dest='mean',type=int, default=None, help='Mean insert. ')
    parser.add_argument('--sd', dest='sd', type=int, default=None, help='Stddev insert.')
    parser.add_argument('--min_size', dest='min_size', type=int, default=None, help='Min insert size (only effective if distr=uniform). ')
    parser.add_argument('--max_size', dest='max_size', type=int, default=None, help='Max insert size (only effective if distr=uniform). ')
   

    args = parser.parse_args()

    if args.distr == 'normal':
        if args.mean == None or args.sd == None:
            print "Argument 'distr' set to 'normal', need both --mean and --sd specified to continue. "
            sys.exit()        
    elif args.distr == 'uniform':
        if args.min_size == None or args.max_size == None:
            print "Argument 'distr' set to 'uniform', need both --min_size and --max_size specified to continue. "
            sys.exit()
    elif args.distr == 'mix':
        if args.mean == None or args.sd == None:
            print "Argument 'distr' set to 'mix', need both --mean and --sd specified to continue. "
            sys.exit()
    else:
        print "Argument 'distr' only takes 'normal' or 'uniform'. You specified: {0}. Exiting without continue. ".format(args.distr)
        sys.exit()
    if (args.contigs and not args.genome) or ( args.genome and not args.contigs):
        print "Either both or none of --genome and --contigs needs to be specified. Exiting without continue. "
        sys.exit()

    main(args)








