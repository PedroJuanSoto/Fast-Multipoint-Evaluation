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
	U = fft(poly_make_deg_eq(p1, d)[::-1])
	V = fft(poly_make_deg_eq(p2, d)[::-1])
	res = list(ifft(U*V).real)
	return [round(x,10) for x in res[::-1]]

#multiply a polynomial by a scalar
def poly_scalar_mult(a, p):
	return [a*pi for pi in p]

#subtract two polynomials
def poly_sub(u, v):
	d = max(len(u), len(v))
	return poly_stand_form([a-b for a,b in zip(poly_make_deg_eq(u, d), poly_make_deg_eq(v, d))])

#add two polynomials
def poly_add(u, v):
	d = max(len(u), len(v))
	return poly_stand_form([a+b for a,b in zip(poly_make_deg_eq(u, d), poly_make_deg_eq(v, d))])

#pad the polynomial p with zeros so that it has new degree = d
def poly_make_deg_eq(p, d):
	return [0] * (d-len(p)) + list(p)

#gets rid of the leading zeroes of p
def poly_stand_form(p):
	for i,a in enumerate(p):
		if a != 0:
			return p[i:]
	if 0 in p: 
		return [0]
	else:
		return []	

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

def newt_inv(p,i):
	g = [1]
	r = 1
	while 2**r < i:
		r += 1 
	for j in range(1,r+1):
		l = poly_scalar_mult(2,g)
		u = poly_mult(g,g)
		m = poly_mult(p,u)
		g = mod_by_x_to_the_i(poly_sub(l,m),2**j)
	return mod_by_x_to_the_i(g, i)

def mod_by_x_to_the_i(p,i):
	return p[-i:]	
	
#print(mod_by_x_to_the_i(poly_mult([1,2,1],newt_inv([1,2,1],3)),3))
#print(rev([[1,2,3],[4,5,6]], 1))
