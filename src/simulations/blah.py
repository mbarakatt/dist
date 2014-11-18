def fibo(seq=[0,1]):
	if seq[-1] > 4000000:
		return seq[0:-1]
	else:
		seq.append(seq[-2] + seq[-1])
		return fibo(seq=seq)


print fibo()

