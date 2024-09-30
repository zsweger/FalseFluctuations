#!/usr/bin/python

import random
import sys

random.seed()
num = int(sys.argv[1])

no_fp = False
try:
    fp = open('rdlist')
except IOError, e:
    no_fp = True

oldnum = []
if not no_fp:
    for line in fp:
        oldnum.append(int(line))

numbers = []
for i in xrange(num):
	while True:
		x = random.randint(0, 2**29)
		if x not in numbers and x not in oldnum:
			numbers.append(x)
			break

output='\n'.join([str(i) for i in numbers])
print(output)
