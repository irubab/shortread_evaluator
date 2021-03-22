import time
import pandas as pd
import ast
import numpy as np
import matplotlib
import sys
import math

matplotlib.use('Agg')

import matplotlib.pyplot as plt

start = time.time()

# read long and short reads
lreads_df = pd.read_csv(snakemake.input[0])

random_sreads = pd.read_csv(snakemake.input[1])

# sample n random short reads (omit this when using all short reads)
n_reads = 1000
random_sreads = random_sreads.sample(n=n_reads)

# Set the limit on x-axis
xlim = 14000

outlier_start_tf = []
outlier_end_tf = []
outlier_overlap_tf = []

normal_start_tf = []
normal_end_tf = []
normal_overlap_tf = []

overlap = {}

start_deviation_list = []
end_deviation_list = []
percentage_list = []

zero_overlap = 0
hundred_overlap = 0
temp_reads_x = []
temp_reads_y = []
temp_reads_z = []

malformed = []


def cigar_parse_start(cigar_tuples, start_position):
    counter = 0
    start_position_cigar = 0
    remainder = 0
    current_tuple_index = 0
    for j in range(len(cigar_tuples)):
        if cigar_tuples[j][0] == 0 or cigar_tuples[j][0] == 2:  # Handle 'H' and 'S'
            counter += cigar_tuples[j][1]
            if counter > start_position:
                if cigar_tuples[j][0] == 0:
                    start_position_cigar += start_position - (counter - cigar_tuples[j][1])
                    # Clipped length of remaining operator
                    remainder = counter - start_position
                    # Store the index of the next cigartuple
                    current_tuple_index = j
                    break
            else:
                if cigar_tuples[j][0] == 0:
                    start_position_cigar += cigar_tuples[j][1]
                elif cigar_tuples[j][0] == 2:
                    continue
        elif cigar_tuples[j][0] == 1:
            start_position_cigar += cigar_tuples[j][1]
            continue
    return start_position_cigar, remainder, current_tuple_index


def cigar_parse_end(cigar_tuples, stop_position, remainder):
    counter = 0
    end_position_cigar = 0
    counter += remainder
    if counter > stop_position:
        end_position_cigar = stop_position
    else:
        for j in range(len(cigar_tuples)):
            if cigar_tuples[j][0] == 0 or cigar_tuples[j][0] == 2:  # Handle 'H' and 'S'
                counter += cigar_tuples[j][1]
                if counter > stop_position:
                    if cigar_tuples[j][0] == 0:
                        end_position_cigar += stop_position - (counter - cigar_tuples[j][1])
                        break
                else:
                    if cigar_tuples[j][0] == 0:
                        end_position_cigar += cigar_tuples[j][1]
                    elif cigar_tuples[j][0] == 2:
                        continue
            elif cigar_tuples[j][0] == 1:
                end_position_cigar += cigar_tuples[j][1]
                continue
    return end_position_cigar

def calculate_percentage_overlap(short_ref_pos, long_ref_pos, short_read_length):
    if len(short_ref_pos) == 1:
        percentage_overlap = 0
    else:
        percentage_overlap = len(list(set(reference_positions) & set(np.arange(ref_pos[0], ref_pos[0] + ref_pos[-1])))) / short_read_length * 100
    return percentage_overlap

reads = 0

for i in random_sreads.index:
    temp = []
    counter = 0
    reference_positions = []
    ref_pos = ast.literal_eval(random_sreads['ref_positions'][i])


    x = random_sreads['qname'][i].split('/')
    lread_name = x[0] + '/' + x[1] + '/ccs'
    start_position_temp = x[2].split('_')
    start_position = start_position_temp[1] if random_sreads['read_type'][i] == 'R1' else start_position_temp[3]

    y = lreads_df.loc[lreads_df['qname'] == lread_name]

    long_index = y.index[0]

    try:
        actual_start = random_sreads['ref_start'][i]
        actual_end = random_sreads['ref_end'][i]
        if math.isnan(actual_start) or math.isnan(actual_end):
            continue
        short_in_long_position = 0
        # read orientation: FF
        if lreads_df['read_alignment'][long_index] == False and random_sreads['read_alignment'][i] == False:
            reads += 1

            cigar_tuples = ast.literal_eval(lreads_df['cigar_tuples'][long_index])
            start_position_cigar, remainder, current_tuple_index = cigar_parse_start(cigar_tuples, int(start_position))
            end_position_cigar = cigar_parse_end(cigar_tuples[current_tuple_index+1:], random_sreads['qlength'][i], remainder)

            expected_start_w_cigar = lreads_df['ref_start'][long_index] + start_position_cigar
            expected_end_w_cigar = expected_start_w_cigar + end_position_cigar
            # print(expected_end_w_cigar, actual_end)

            start_deviation = abs(expected_start_w_cigar - actual_start)
            end_deviation = abs(expected_end_w_cigar - actual_end)

            reference_positions.extend(range(expected_start_w_cigar, expected_end_w_cigar))
            percentage_overlap = calculate_percentage_overlap(ref_pos, reference_positions, random_sreads['qlength'][i])

        # read orientation: FR
        elif lreads_df['read_alignment'][long_index] == False and random_sreads['read_alignment'][i] == True:
            reads+=1

            cigar_tuples = ast.literal_eval(lreads_df['cigar_tuples'][long_index])
            start_position_cigar, remainder, current_tuple_index = cigar_parse_start(cigar_tuples, int(start_position))
            end_position_cigar = cigar_parse_end(cigar_tuples[current_tuple_index + 1:], random_sreads['qlength'][i], remainder)

            expected_start_w_cigar = lreads_df['ref_start'][long_index] + start_position_cigar - random_sreads['qlength'][i]
            expected_end_w_cigar = lreads_df['ref_start'][long_index] + start_position_cigar

            start_deviation = abs(expected_start_w_cigar - actual_start)
            end_deviation = abs(expected_end_w_cigar - actual_end)


            reference_positions.extend(range(expected_start_w_cigar, expected_end_w_cigar))
            percentage_overlap = calculate_percentage_overlap(ref_pos, reference_positions, random_sreads['qlength'][i])

        # read orientation: RF
        elif lreads_df['read_alignment'][long_index] == True and random_sreads['read_alignment'][i] == False:
            reads+=1

            cigar_tuples = ast.literal_eval(lreads_df['cigar_tuples'][long_index])[::-1]
            start_position_cigar, remainder, current_tuple_index = cigar_parse_start(cigar_tuples, int(start_position))
            end_position_cigar = cigar_parse_end(cigar_tuples[current_tuple_index + 1:], random_sreads['qlength'][i], remainder)

            expected_start_w_cigar = lreads_df['ref_end'][long_index] - start_position_cigar
            expected_end_w_cigar = expected_start_w_cigar + random_sreads['qlength'][i]
            # print(expected_start_w_cigar, expected_end_w_cigar, actual_start, actual_end)

            start_deviation = abs(expected_start_w_cigar - actual_start)
            end_deviation = abs(expected_end_w_cigar - actual_end)

            reference_positions.extend(range(int(expected_start_w_cigar), int(expected_end_w_cigar)))
            percentage_overlap = calculate_percentage_overlap(ref_pos, reference_positions, random_sreads['qlength'][i])


        # read orientation: RR
        elif lreads_df['read_alignment'][long_index] == True and random_sreads['read_alignment'][i] == True:
            reads+=1

            cigar_tuples = ast.literal_eval(lreads_df['cigar_tuples'][long_index])[::-1]
            start_position_cigar, remainder, current_tuple_index = cigar_parse_start(cigar_tuples, int(start_position))
            end_position_cigar = cigar_parse_end(cigar_tuples[current_tuple_index + 1:], random_sreads['qlength'][i], remainder)

            expected_start_w_cigar = lreads_df['ref_end'][long_index] - start_position_cigar - end_position_cigar
            expected_end_w_cigar = lreads_df['ref_end'][long_index] - start_position_cigar
            # print(expected_start_w_cigar, expected_end_w_cigar, actual_start, actual_end)

            start_deviation = abs(expected_start_w_cigar - actual_start)
            end_deviation = abs(expected_end_w_cigar - actual_end)

            reference_positions.extend(range(int(expected_start_w_cigar), int(expected_end_w_cigar)))
            percentage_overlap = calculate_percentage_overlap(ref_pos, reference_positions, random_sreads['qlength'][i])

        else:
            continue

        ########################################## Plotting ########################################

        if start_deviation > xlim:
            # temp.append(110 + 1)  # Set deviation to a number which is 1 greater than the short read length
            start_deviation_list.append(xlim + 1)
            outlier_start_tf.append(random_sreads['qname'][i])
        else:
            start_deviation_list.append(start_deviation)
            normal_start_tf.append(random_sreads['qname'][i])

        if end_deviation > xlim:
            # temp.append(110 + 1)  # Set deviation to a number which is 1 greater than the short read length
            end_deviation_list.append(xlim + 1)
            outlier_end_tf.append(random_sreads['qname'][i])
        else:
            end_deviation_list.append(end_deviation)
            normal_end_tf.append(random_sreads['qname'][i])

        if percentage_overlap > 100:
            percentage_list.append(0)  # No overlap
            outlier_overlap_tf.append(random_sreads['read_alignment'][i])
        else:
            percentage_list.append(percentage_overlap)
            normal_overlap_tf.append(random_sreads['read_alignment'][i])

        overlap[random_sreads['qname'][i]] = [start_deviation, end_deviation, percentage_overlap]
        if percentage_overlap == 100:
            temp_reads_x.append(random_sreads['qname'][i])
            hundred_overlap+=1
        elif percentage_overlap == 0:
            temp_reads_y.append(random_sreads['qname'][i])
            zero_overlap+=1
        elif percentage_overlap < 100 and percentage_overlap > 0:
            temp_reads_z.append(random_sreads['qname'][i])
        

    except ValueError:
        malformed.append(random_sreads['qname'][i])
        continue

print(len(malformed))
print(reads)
print('100% overlap:',(hundred_overlap/n_reads)*100)
print('Zero overlap:',(zero_overlap/n_reads)*100)


fig = plt.figure(figsize=(15,15))

ax1 = fig.add_subplot(3,1,1)
ax1.hist(start_deviation_list, bins = 'fd', label='Start deviation')
ax1.set_title('Start deviation')
ax1.set_xlabel('No. of reads')
ax1.set_ylabel('count of')

ax2 = fig.add_subplot(3,1,2)
ax2.hist(end_deviation_list, bins = 'fd', label='End deviation')
ax2.set_title('End deviation')
ax2.set_xlabel('No. of reads')
ax2.set_ylabel('count of')

ax3 = fig.add_subplot(3,1,3)
ax3.hist(percentage_list, bins = 'fd', label='Percentage Overlap')
ax3.set_title('Percentage Overlap')
ax3.set_xlabel('No. of reads')
ax3.set_ylabel('percentage')

plt.grid(True)
plt.tight_layout()
plt.savefig(snakemake.output[0])

end = time.time()

print('It took: ' + str(end - start) + ' secs')
