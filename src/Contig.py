'''
Created on Sep 29, 2011

@author: ksahlin
'''
class contig:
    def __init__(self, contig_name,contig_scaffold= None ,contig_direction = None ,contig_position = None, contig_length = None,contig_sequence=None, contig_links={} ):
        self.name = contig_name
        self.scaffold = contig_scaffold
        self.direction = contig_direction
        self.position = contig_position
        self.length = contig_length
        self.sequence = contig_sequence
        self.links = contig_links
        #,contig_leftLinks=None,contig_rightLinks=None 
        
#    def UpdateIndex(self,contig_position,contig_direction,contig_scaffold):
#        return()