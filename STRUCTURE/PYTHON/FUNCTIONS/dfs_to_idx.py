# Python script for converting PLUMED output to idx file
# Optional command line args: (1) dfs_file [dfs_surf.dat], (2) idx_file [idx.dat]

import sys

dfs_file = 'dfs_surf.dat'
idx_file = 'idx.dat'

if len(sys.argv) > 1:
    dfs_file = sys.argv[1]
    if len(sys.argv) > 2:
        idx_file = sys.argv[2]

with open(dfs_file, 'r') as f: lines = f.readlines()
with open(idx_file, 'w+') as o:
    for i in range(len(lines)):
        if i % 2:
            for idx in lines[i][18:].split(): o.write('{} '.format(int(idx)+1))
            o.write('\n')
        else: o.write('{} '.format(lines[i].split()[-1]))