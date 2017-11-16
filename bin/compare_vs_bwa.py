#!/usr/bin/python3

import sys
# Requires 2 arguments
if len(sys.argv) != 3:
    print("Usage: ./compare_vs_bwa.py <our mapper output> <bwa output>")
    print("Requires two arguments. Arguments provided: {}".format(len(sys.argv)))
    sys.exit()

with open(sys.argv[1]) as f:
    our_output = f.read().split("\n")

with open(sys.argv[2]) as f:
    bwa_output = f.read().split("\n")

reverse = 16
unmapped = 4

false_positive = 0
false_negative = 0
total_mappings = 0
total_other = 0
sequence = ""
normalize = 5
for i in range(len(reference)):
    # BWA: first item is the name
    # Other: i % 3 == 0 is the name
    #        i % 3 == 1 is the forward
    #        i % 3 == 2 is the reverse
    if i%2 == 1:
        d1 = [int(s) for s in reference[i].split(' ')[:-1]]
        d2 = [int(s) for s in data2[i].split(' ')[:-1]]
        total_mappings += len(d1)
        total_other += len(d2)
    if i%4 == 0:
        sequence = reference[i]
    if reference[i] == data2[i]:
        continue

    d1 = [int(int(s)/normalize) for s in reference[i].split(' ')[:-1]]
    d2 = [int(int(s)/normalize) for s in data2[i].split(' ')[:-1]]
    for d in d1:
        if d not in d2:
            if d+1 not in d2:
                if d-1 not in d2:
                    false_negative += 1
                    #print(sequence)
                    #if i%4 == 3:
                    #    print("reverse")
                    #print(d*normalize)
    for d in d2:
        if d not in d1:
            if d+1 not in d1:
                if d-1 not in d1:
                    false_positive += 1

print('Errors: {}p {}n = {}. Total Mappings = {}. Total other = {}'.format(false_positive, false_negative, false_positive + false_negative, total_mappings, total_other))
print('False Negative Rate: {}'.format(0 if false_negative == 0 else false_negative/total_mappings))
print('False Positive Rate: {}'.format(0 if false_positive == 0 else false_positive/total_mappings))
