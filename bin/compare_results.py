#!/usr/bin/python3

import sys
# Requires 2 arguments
if len(sys.argv) != 3:
    print("Requires two arguments. Arguments provided: {}".format(len(sys.argv)))

with open(sys.argv[1]) as f:
    reference = f.read().split("\n")

with open(sys.argv[2]) as f:
    data2 = f.read().split("\n")

if len(reference ) != len(data2):
    print("Lengths do not match.")

false_positive = 0
false_negative = 0
total_mappings = 0
for i in range(len(reference)):
    if i%2 == 1:
        d1 = [int(s) for s in reference[i].split(' ')[:-1]]
        d2 = [int(s) for s in data2[i].split(' ')[:-1]]
        total_mappings += len(d1)
    if reference[i] == data2[i]:
        continue

    d1 = [int(int(s)/100) for s in reference[i].split(' ')[:-1]]
    d2 = [int(int(s)/100) for s in data2[i].split(' ')[:-1]]
    for d in d1:
        if d not in d2:
            false_negative += 1
    for d in d2:
        if d not in d1:
            false_positive += 1

print('Errors: {}p {}n = {}. Total Mappings = {}'.format(false_positive, false_negative, false_positive + false_negative, total_mappings))
print('False Negative Rate: {}'.format(0 if false_negative == 0 else false_negative/total_mappings))
print('False Positive Rate: {}'.format(0 if false_positive == 0 else false_positive/total_mappings))
