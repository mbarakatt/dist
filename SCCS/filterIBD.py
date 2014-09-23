import sys


for line in sys.stdin:
	sp=line.split()
	if "GWAS" in sp[1] and "GWAS" in sp[3]:
		sys.stdout.write(line)

