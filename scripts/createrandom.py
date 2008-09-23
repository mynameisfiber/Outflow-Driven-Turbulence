#!/usr/bin/python

import sys, random

file(sys.argv[1],"w+").write("\n".join(["%.10f"%random.random() for i in range(int(sys.argv[2]))]))