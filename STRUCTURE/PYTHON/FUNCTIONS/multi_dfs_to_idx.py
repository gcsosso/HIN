# Python script for converting PLUMED output to idx file
# Optional command line args: (1) dfs_file prefix [cls], (2) number of dat files [100],
#                             (3) idx_file [idx.dat], (4) min cluster size for inclusion [40]
# Will read cls_1.dat, cls_2.dat, ... , cls_100.dat

import sys

dfs_file = 'cls'
N = 100
idx_file = 'idx.dat'
min_size = 40

if len(sys.argv) > 1:
    dfs_file = sys.argv[1]
    if len(sys.argv) > 2:
        N = int(sys.argv[2])
        if len(sys.argv) > 3:
            idx_file = sys.argv[3]
            if len(sys.argv) > 4:
                min_size = int(sys.argv[4])

with open('{}_{}.dat'.format(dfs_file,N), 'r') as f: lines = f.readlines()
lentraj = int(len(lines)/2)
with open(idx_file, 'w+') as o:
    for i in range(lentraj):
        line_count = 0
        indices = ''
        for j in range(1,N+1):
            with open('{}_{}.dat'.format(dfs_file,j), 'r') as f: lines = f.readlines()
            if int(lines[2*i].split()[-1]) >= min_size:
                line_count += int(lines[2*i].split()[-1])
                for idx in lines[2*i+1][18:].split(): indices = '{} {}'.format(indices, int(idx)+1)
            
        o.write('{}{}\n'.format(line_count,indices))
