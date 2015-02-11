'''
Created on Sep 23, 2011

@author: ksahlin
'''

import pysam
import Contig,Scaffold
import networkx as nx
from collections import defaultdict
import os,sys

##
# Opens a .bam or .sam file and returns the file
# object.
#
# @param bam_file_path Path to the .bam or .sam.
# 
# @return File object for .bam or .sam.
#
def open_bam_file(bam_file_path):
    bam_file_name, bam_file_ext = os.path.splitext(bam_file_path)
    if bam_file_ext == ".bam":
        return pysam.Samfile(bam_file_path, 'rb')
    elif bam_file_ext == ".sam":
        return pysam.Samfile(bam_file_path, 'r')
    else:
        return IOError("open_bam_file: File must be either .bam or .sam.")

def is_proper_aligned_unique_innie(read):
    return (read.is_reverse and not read.mate_is_reverse and  read.tlen < 0 and read.rname == read.mrnm) or \
                (not read.is_reverse and read.mate_is_reverse and read.is_read2 and read.tlen > 0 and read.rname == read.mrnm ) \
                and not read.mate_is_unmapped and not read.is_unmapped and read.mapq > 10 and not read.is_secondary
def is_proper_aligned_unique_outie(read):
    return (read.is_reverse and not read.mate_is_reverse and  read.tlen > 0 and read.rname == read.mrnm) or \
                (not read.is_reverse and read.mate_is_reverse and read.is_read2 and read.tlen < 0 and read.rname == read.mrnm ) \
                and not read.mate_is_unmapped and not read.is_unmapped and read.mapq > 10 and not read.is_secondary
def is_unique_read_link(read):
    # if  not read.is_unmapped and not read.mate_is_unmapped and read.rname != read.mrnm \
    # and read.opt('XT')=='U' and not read.is_secondary and read.rlen != read.alen:
    #     print read
    return not read.is_unmapped and not read.is_secondary and read.mapq >= 10



class BamParser(object):
    """docstring for BamParser"""
    def __init__(self, bam_file):
        super(BamParser, self).__init__()
        self.bam_file = open_bam_file(bam_file)

        self.contig_lengths = dict(zip(self.bam_file.references,self.bam_file.lengths))   

    def aligned_reads(self,aligner):
        if aligner == 'bwa' or aligner == 'bwa_mem':
            for read in self.bam_file:
                if not read.is_unmapped:
                    yield read 
            self.bam_file.seek(0)
        elif aligner == 'bowtie':
            for read1,read2 in zip(self.bam_file,self.bam_file2):
                if not read1.is_unmapped and not read2.is_unmapped:
                    yield read1,read2           
                elif not read1.is_unmapped:
                    yield read1,False
                elif not read2.is_unmapped:
                    yield False, read2

    def unique_reads_on_different_references(self):
        read_pairs = {}
        for read in self.bam_file:
            if is_unique_read_link(read):
                #tmp for tests:
                #print read.qname[:-1]
                #print read_pairs
                if read.qname in read_pairs:
                    #print 'lol'
                    read2 = read_pairs[read.qname]
                    if read.tid != read2.tid:    
                        yield read, read_pairs[read.qname]
                    #else: 
                    #    pass
                    del read_pairs[read.qname]
                else:
                    read_pairs[read.qname] = read

        self.bam_file.seek(0)


def PE(Contigs,Scaffolds,bamfile,mean,std_dev,scaffold_indexer,F,read_len):
    G=nx.Graph()
    #print 'Parsing BAM file...'
    #read_len=50
    #informative_pair={81:(False,True),97:(True,False),113:(False,False),65:(True,True)}
    #I switched to look at mates instead since BWA can give false flag combinations for
    # read-mate when read is mapped but not mate eg 97-149 81-165. But the reverse
    #does not happen.
    #informative_pair={161:(True,False),145:(False,True),129:(True,True),177:(False,False)}
    #threshold=800
    def AddEdges(Contigs,Scaffolds,bamfile,mean,std_dev,scaffold_indexer,F,read_len):
        #Clean contig_library
        bam_object = BamParser(bamfile)
        singeled_out=0
        cont_lengths= bam_object.bam_file.lengths
        cont_lengths=[int(nr) for nr in cont_lengths]  #convert long to int object
        #print cont_lengths
        cont_names = bam_object.bam_file.references
        ####### WHEN ADDING SHORTER CONTIGS NOT INCLUDED IN THE SCAFFOLDING, 
        ####### WE NEED TO ALSO INITIALIZE OBJECTS FOR THESE, THIS SHOULD BE DONE SOMEWHERE HERE
        for i in range(0,len(cont_names)):
            if cont_lengths[i] >= 300:
                C=Contig.contig(cont_names[i])   # Create object contig
                C.length = cont_lengths[i]
                C.scaf_length = C.length        # Initially, scaffold consists of only this contig
                C.direction = True              # always in same direction first, False=reverse
                C.position = 0                  #position always 0
                C.links = {}
                Contigs[C.name] = C              # Create a dict with name as key and the object container as value
                S=Scaffold.scaffold('s'+str(scaffold_indexer),[C],C.length)  # Create object scaffold
                Scaffolds[S.name]=S
                C.scaffold=S.name
                G.add_node((S.name,'L'),length=cont_lengths[i])
                G.add_node((S.name,'R'),length=cont_lengths[i])
                scaffold_indexer+=1
        
        
        #Create "node graph" of contigs (that passed the length criteria). Each having a left and right node
        #print 'Nr of contigs/scaffolds included in scaffolding: '+ str(len(Scaffolds))#,Scaffolds.keys()
        
        for scaffold_ in Scaffolds:
            G.add_edge((scaffold_,'L'),(scaffold_,'R'),nr_links=None)    #this is a scaffold object but can be both a single contig or a scaffold.
        
        
        # Create the link edges in the graph by fetching info from bam file


        def nr_softclipps(read):
            max_soft = 0
            for type_,length in read.cigar:
                if type_ == 4 and length >= max_soft:
                    max_soft = length
            return max_soft

        global_max_softclipps = 0
        global_min_obs = 100000 
        links_used = 0
        #r_len = float(read_len)
        for read1,read2 in bam_object.unique_reads_on_different_references():
            contig1=bam_object.bam_file.getrname(read1.rname)
            contig2=bam_object.bam_file.getrname(read2.rname)
            max_soft_readpair = max(nr_softclipps(read1),nr_softclipps(read2))
            if max_soft_readpair > global_max_softclipps:
                global_max_softclipps = max_soft_readpair
            # print read1.cigar
            #if read1.qlen/r_len < 0.7 or read2.qlen/r_len < 0.7:
            #    continue
            #     print 'midddle1',o1, o1+o2, read1.pos, read1.mapq,read1.qlen,read1.rlen, read1.cigar, read1.tags
            # if read2.qlen < 50:
            #     print 'midddle2',o2, o1+o2, read2.pos, read2.mapq, read2.qlen,read2.rlen, read2.cigar, read2.tags
            if contig1 in Contigs and contig2 in Contigs:                
                (read_dir,mate_dir) = (not read1.is_reverse,not read2.is_reverse )
                scaf1=Contigs[contig1].scaffold
                scaf2=Contigs[contig2].scaffold                    
                #Calculate actual position on scaffold here
                #position1 cont/scaf1
                cont_dir1 = Contigs[contig1].direction  #if pos : L if neg: R
                cont1_pos = Contigs[contig1].position
                readpos = read1.pos
                cont1_len = Contigs[contig1].length
                s1len = Scaffolds[scaf1].s_length
                #position1 cont1/scaf1                        
                cont_dir2 = Contigs[contig2].direction
                cont2_pos = Contigs[contig2].position
                matepos = read2.pos
                cont2_len = Contigs[contig2].length
                s2len = Scaffolds[scaf2].s_length 
                (obs,scaf_side1,scaf_side2, (o1,o2))=PosDirCalculatorPE(cont_dir1,read_dir,cont1_pos,readpos,s1len,cont1_len,cont_dir2,mate_dir,cont2_pos,matepos,s2len,cont2_len,read_len) 
                if obs < mean+ 4*std_dev: 
                    links_used += 1
                    if (scaf2,scaf_side2) not in G[(scaf1,scaf_side1)]:
                        G.add_edge((scaf2,scaf_side2),(scaf1,scaf_side1),nr_links=1,gap_dist=[obs],obs_pos=set() )
                        G[(scaf2,scaf_side2)][(scaf1,scaf_side1)]['obs_pos'].add((o1,o2))
                        if o1 < global_min_obs:
                            global_min_obs = o1
                        if o2 < global_min_obs:
                            global_min_obs = o2 
                    #print 'Added edge'
                    else:
                        try:
                            if (o1,o2) in G.edge[(scaf1,scaf_side1)][(scaf2,scaf_side2)]['obs_pos']:
                                continue
                        except KeyError:
                            #print G.edge[(scaf1,scaf_side1)][(scaf2,scaf_side2)]
                            continue

                        # if (o1,o2) in G.edge[(scaf1,scaf_side1)][(scaf2,scaf_side2)]['obs_pos']:
                        #     #print 'detected duplicate'
                        #     continue
                        else:
                            G.edge[(scaf1,scaf_side1)][(scaf2,scaf_side2)]['nr_links'] += 1
                            G.edge[(scaf1,scaf_side1)][(scaf2,scaf_side2)]['gap_dist'].append(obs)  
                            G.edge[(scaf1,scaf_side1)][(scaf2,scaf_side2)]['obs_pos'].add((o1,o2))  
                            G.edge[(scaf1,scaf_side1)][(scaf2,scaf_side2)]['obs_pos'].add((o2,o1))  

                            if o1 < global_min_obs:
                                global_min_obs = o1
                            if o2 < global_min_obs:
                                global_min_obs = o2
                            # if o1 < 50:
                            #     print o1, o1+o2, read1.pos, read1.mapq,read1.qlen,read1.rlen, read1.cigar, read1.tags
                            #     #print fancy_str(read1)
                            # if o2 < 50:
                            #     print o2, o1+o2, read2.pos, read2.mapq, read2.qlen,read2.rlen, read2.cigar, read2.tags
                            #     #print fancy_str(read2)                                
                                

        print 'Max softclipps:', global_max_softclipps
        print 'Min obs:', global_min_obs
        # sys.exit()
        #print 'Nr links used:', links_used
        return global_max_softclipps
    
    max_softclipps = AddEdges(Contigs,Scaffolds,bamfile,mean,std_dev,scaffold_indexer,F,read_len)

    return(G,Contigs,Scaffolds,F,scaffold_indexer,max_softclipps)

def RemoveBugEdges(G,fishy_edges):
    edges_removed = 0
    for edge_tuple,nr_links in fishy_edges.items():             
        if edge_tuple[1] in G and edge_tuple[0] in G[edge_tuple[1]]:
            if nr_links >= G[edge_tuple[0]][edge_tuple[1]]['nr_links']:
                G.remove_edge(edge_tuple[0],edge_tuple[1])  
                edges_removed += 1 
#print 'Number of BWA buggy edges removed: ', edges_removed           
    return()

def CheckDir(cont_obj1,cont_obj2,alignedread):
    (read_dir,mate_dir) = (not alignedread.is_reverse,not alignedread.mate_is_reverse )
    cont_dir1 = cont_obj1.direction  #if pos : L if neg: R
    #position2 cont2/scaf2                        
    cont_dir2 = cont_obj2.direction
    (gap,scaf_side1,scaf_side2) = PosDirCalculatorPE(cont_dir1,read_dir,0,0,0,0,cont_dir2,mate_dir,0,0,0,0,0)    
    return(scaf_side1,scaf_side2)

def PosDirCalculatorPE(cont_dir1,read_dir,cont1pos,readpos,s1len,cont1_len,cont_dir2,mate_dir,cont2pos,matepos,s2len,cont2_len,read_len):
    #calculates the positions of the two reads on theis respective contig.
    #o_1^i is denoted gap1 in the code (a bit misleading)
    #o_2^i is denoted gap2 in the code (a bit misleading)
    # o^i is then o_1^i+o_2^i and is denoted gap here
    if cont_dir1 and read_dir:
        gap1=s1len-cont1pos-readpos
        read_side1='R'
    if cont_dir2 and mate_dir:
        gap2=s2len-cont2pos-matepos
        read_side2='R'
    if (not cont_dir1) and read_dir:
        gap1=cont1pos+(cont1_len-readpos)
        read_side1='L'
    if (not cont_dir2) and mate_dir:
        gap2=cont2pos+(cont2_len-matepos)
        read_side2='L'
    if cont_dir1 and not read_dir:
        gap1=cont1pos + readpos + read_len
        read_side1='L'
    if cont_dir2 and not mate_dir:
        gap2=cont2pos + matepos + read_len
        read_side2='L'
    if not cont_dir1 and not read_dir:
        gap1= s1len - cont1pos - (cont1_len-readpos -read_len)
        read_side1='R'
    if not cont_dir2 and not mate_dir:
        gap2= s2len - cont2pos - (cont2_len-matepos -read_len)
        read_side2='R'
    obs=gap1+gap2
    if read_side1 == 'L':
        scaf_side1 = 'L'
    if read_side2 == 'L':
        scaf_side2 = 'L'
    if read_side1 == 'R':
        scaf_side1 = 'R'
    if read_side2 == 'R':
        scaf_side2 = 'R'
    return(obs,scaf_side1,scaf_side2,(gap1,gap2))







