#Program to process regression output and extract the best models

import sys
import csv
import pandas as pd

data = sys.argv[1]
pass_models = sys.argv[2]

p1 = []
with open(data) as output:
	for line in output.readlines():
		try:
			line = line.strip()
			contents = line.split('\t')
			roc = float(contents[-2])
			f1 = float(contents[-1])

			if(roc>=0.5 and f1>=0.5):
				p1.append(contents)
		except:
			print(line)
			continue

with open(pass_models, 'w') as f1:
	for model in p1:
		print('\t'.join(model), file=f1)

print("Output processed.")
