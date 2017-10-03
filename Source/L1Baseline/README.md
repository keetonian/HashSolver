== Baseline Hash Implementation ==

Converts bases into hash keys, where each base is 2 bits.

A: 00
C: 01
G: 10
T: 11

This hash table uses a 1:1 mapping; the hash is used as an array index.
