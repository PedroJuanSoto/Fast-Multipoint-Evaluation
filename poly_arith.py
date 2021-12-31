from math import ceil, log
from numpy.fft import fft, ifft
import numpy as np
from copy import deepcopy

#FFT polynomial multiplication
def poly_mult(p1, p2):
	if not p1 or not p2:
		return []
	deg1 = len(p1) - 1 
	deg2 = len(p2) - 1
	d = deg1 + deg2 + 1
	U = fft(poly_extend(p1, d)[::-1])
	V = fft(poly_extend(p2, d)[::-1])
	res = list(ifft(U*V).real)
	return [round(x,10) for x in res[::-1]]
	
def poly_extend(p, d):
	return [0] * (d-len(p)) + list(p)

def poly_div(p, q):
	#return np.poly1d(poly_divmod(list(p), list(q))[1])
	return np.polydiv(p, q)[1]

def rev(p, k):
	l, r = p
	rev_p = [deepcopy(r),deepcopy(l)] 
	rev_p[0].reverse()
	rev_p[1].reverse()
	for i in range(k):
		rev_p[0].append(rev_p[1].pop(0))
	return rev_p


#print(rev([[1,2,3],[4,5,6]], 1))
