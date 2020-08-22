#!/usr/bin/env python

import sys

in_file = open(sys.argv[1], "r")

out_file = open(sys.argv[2], "w")

count = 0
for line in in_file:
        out_file.write(line.replace("@GClone|", ("@GClone_"+str(count)+"|")))
        count += line.count("@GClone|")

in_file.close()
out_file.close()
