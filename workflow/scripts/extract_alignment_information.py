import time
import pandas as pd
import pysam
import sys

# Calculate the % overlap between alignments of short and long reads w.r.t a reference genome.

start = time.time()

short_samfile = pysam.AlignmentFile(snakemake.input[0], 'r' ) # Access these files from the server
long_samfile = pysam.AlignmentFile(snakemake.input[1],'r')  # Access these files from the server


long_data = {'read_alignment':[], 'qname':[],'query':[],'qlength':[],'qalign_start':[],'qalign_end':[],'ref_id':[],'ref_name':[],'ref_start':[],'ref_end':[],'map_qual':[], 'cigar_tuples':[]}

short_data = {'read_alignment':[], 'qname':[], 'query':[],'qlength':[],'qalign_start':[],'qalign_end':[],'ref_id':[],'ref_name':[],'ref_start':[],'ref_end':[],'map_qual':[], 'cigar_tuples':[], 'ref_positions':[], 'read_type':[]}


for l in long_samfile.fetch():
    long_data['read_alignment'].append(l.is_reverse)
    long_data['qname'].append(l.query_name)
    long_data['query'].append(l.query_sequence)
    long_data['qlength'].append(l.query_length)
    long_data['qalign_start'].append(l.query_alignment_start)
    long_data['qalign_end'].append(l.query_alignment_end)
    long_data['ref_id'].append(l.reference_id)
    long_data['ref_name'].append(l.reference_name)
    long_data['ref_start'].append(l.reference_start)
    long_data['ref_end'].append(l.reference_end)
    long_data['map_qual'].append(l.mapping_quality)
    long_data['cigar_tuples'].append(l.cigartuples)

lreads_df = pd.DataFrame.from_dict(long_data)
lreads_df.to_csv(snakemake.output[0], index=False) 
long_samfile.close()


for s in short_samfile.fetch():
    short_data['read_alignment'].append(s.is_reverse)
    if s.is_read1:
        short_data['read_type'].append('R1')
    elif s.is_read2:
        short_data['read_type'].append('R2')
    short_data['qname'].append(s.query_name)
    short_data['query'].append(s.query_sequence)
    short_data['qlength'].append(s.query_length)
    short_data['qalign_start'].append(s.query_alignment_start)
    short_data['qalign_end'].append(s.query_alignment_end)
    short_data['ref_id'].append(s.reference_id)
    short_data['ref_name'].append(s.reference_name)
    short_data['ref_start'].append(s.reference_start)
    short_data['ref_end'].append(s.reference_end)
    short_data['map_qual'].append(s.mapping_quality)
    short_data['cigar_tuples'].append(s.cigartuples)
    ref_align = len(s.get_reference_positions())
    temp = []
    if ref_align == 0:
        temp = [0,0,0]
    else:
        temp.append(s.get_reference_positions()[0])
        temp.append(s.get_reference_positions()[-1])
        temp.append(ref_align)
    short_data['ref_positions'].append(temp)

sreads_df = pd.DataFrame.from_dict(short_data)
sreads_df.to_csv(snakemake.output[1], index=False)
short_samfile.close()

end = time.time()
