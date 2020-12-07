import sys
from collections import defaultdict
gene_fdr = defaultdict(list)
gene_fc = defaultdict(list)
with open(sys.argv[1]) as f:
	for line in f:
		if 'name' in line:
			continue
		else:
			line = line.strip().split('\t')
			gene = str(line[3])
			fdr = float(line[15])
			fc = float(line[17])
		gene_fdr[gene].append(fdr)
		gene_fc[gene].append(fc)

for k,v in gene_fdr.items():
	print(k,sum(v)/len(v),sum(gene_fc[k])/len(gene_fc[k]))
