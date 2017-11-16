#!/usr/bin/python3
import sys

if len(sys.argv) != 3:
    print("Usage: annotate.py <data structures> <trace file>")
    sys.exit()

with open(sys.argv[1]) as f:
    data = f.read().split('\n')

ignore = ["clock", "decompress", "convert2bit1"]

if data[-1] == "":
    data = data[:-1]
names = [s.split(":")[0] for s in data if s]
beginning = [int(s.split(":")[1].split("-")[0], 16) for s in data if s]
ending = [int(s.split(":")[1].split("-")[1], 16) for s in data if s]

# Read Traces
with open(sys.argv[2]) as f:
    data = f.read().split('\n')

total_ops = 0
skip = False
for d in data:
    # Function Names
    if len(d) == 0:
        continue
    if d[0] != " ":
        skip = False
        for i in ignore:
            if i in d:
                skip = True
        if skip:
            continue
        print("\t{}".format(total_ops))
        total_ops = 0
        print(d)
        continue
    if skip:
        continue
    
    # Memory Traces
    data_parts = d.split("\t")
    address = int(data_parts[-2], 16)
    access_type = data_parts[-3]
    num_ops_before = int(data_parts[-4])
    for i in range(len(names)):
        if address >= beginning[i] and address <= ending[i]:
            print("\t{}\t{}\t{} {}".format(total_ops, access_type, address, names[i]))
            total_ops = 0
            break
    if(i == len(names)-1):
        total_ops += num_ops_before

