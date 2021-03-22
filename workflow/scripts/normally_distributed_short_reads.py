import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

import numpy as np
from scipy.stats import uniform, norm
from scipy.stats import truncnorm
import re

forward_sequences = []
reverse_sequences = []
read_counter = 0

all_sequences = list(SeqIO.parse(snakemake.input[0], 'fastq'))


def get_truncated_norm(m=0, sd=1, low=0, upp=10):
    return truncnorm((low - m) / sd, (upp - m) / sd, loc=m, scale=sd)


def sample_short_reads(mean, standard_dev, sreadLen):
    temp = 0
    for record in all_sequences:
        longReadLength = len(record.seq)
        counter = 0
        readCounter = 0
        qual_score = record.letter_annotations['phred_quality']

        while counter <= 0:
            insertSize = get_truncated_norm(m=mean, sd=standard_dev, low=200, upp=300).rvs()
            if (readCounter + 2 * sreadLen + int(insertSize)) <= longReadLength:
                r1 = record.seq[readCounter:readCounter + sreadLen]  # The 0 is used to index the element stored as [[]] in r.
                q1 = qual_score[readCounter:readCounter + sreadLen]

                r2_temp = record.seq[readCounter + sreadLen + int(insertSize): readCounter + sreadLen + int(insertSize) + sreadLen]
                r2 = r2_temp.reverse_complement()
                q2_temp = qual_score[readCounter + sreadLen + int(insertSize): readCounter + sreadLen + int(insertSize) + sreadLen]
                q2 = q2_temp[::-1]

                idr = record.id + '_' + str(readCounter) + '_' + str(insertSize) + '_' + str(readCounter + sreadLen + int(insertSize) + sreadLen)

                forward_sequences.append(SeqRecord(r1, id=idr, name=idr, description=idr, letter_annotations={'phred_quality': q1}))

                reverse_sequences.append(SeqRecord(r2, id=idr, name=idr, description=idr, letter_annotations={'phred_quality': q2}))
                readCounter += sreadLen
                continue
            else:
                counter = 1
                continue

        temp += 1
    return forward_sequences, reverse_sequences


fs, rs = sample_short_reads(250, 1, 110)

# Store the contents into a fastq file

with open(snakemake.output[1], 'w') as fastqfile:
    SeqIO.write(fs, fastqfile, 'fastq-illumina')
with open(snakemake.output[2], 'w') as fastqfile_rev:
    SeqIO.write(rs, fastqfile_rev, 'fastq-illumina')
