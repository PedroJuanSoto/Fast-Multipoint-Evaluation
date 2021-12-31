from math import ceil, log
from numpy.fft import fft, ifft
import numpy as np

#FFT polynomial multiplication
def poly_mult(p1, p2):
	if not p1 or not p2:
		return []
	deg1 = len(p1) - 1 
	deg2 = len(p2) - 1
	d = deg1 + deg2 + 1
	def poly_extend(p, d):
		return [0] * (d-len(p)) + list(p)
	U = fft(poly_extend(p1, d)[::-1])
	V = fft(poly_extend(p2, d)[::-1])
	res = list(ifft(U*V).real)
	return [round(x,10) for x in res[::-1]]
	
#Going up the product tree in the multipoint evaluation algorithm 
def create_tree(r):
	num = len(r)
	tree = [[]]
	for i in range(num):
		tree[0].append([1,-r[i]]) 
	i = 1
	j = 2
	while j <= num:
		tree.append([])
		for k in range(num // j):
	  		prod = poly_mult(tree[i-1][k*2], tree[i-1][k*2 + 1])    
	  		tree[i].append(prod) 
		j = j*2
		i = i+1
	return tree 
	
def poly_div(p, q):
	#return np.poly1d(poly_divmod(list(p), list(q))[1])
	return np.polydiv(p, q)[1]
		
#Going down the product tree 
def mult_eval(poly, points):
	num = len(points)
	if num == 1:
		return poly
	j = 1
	while j < num:
		j = 2*j
	eval_points = []
	if num == j:
		for i in range(num):
	  		eval_points.append(points[i])
	else:
		for i in range(num):
	  		eval_points.append(points[i])
		for i in range(j - num):
	  		eval_points.append(0)
	tree = create_tree(eval_points)
	r0 = poly_div(poly, tree[len(tree)-2][0])
	r1 = poly_div(poly, tree[len(tree)-2][1])
	list1 = mult_eval_recurse(r0, eval_points[:j//2], tree, 1, 0)
	list2 = mult_eval_recurse(r1, eval_points[j//2:], tree, 1, 1)
	return (list1 + list2)[:num] 

#A helper recursion function for going down product tree
def mult_eval_recurse(poly, points, tree, depth, index):
	num = len(points)
	if num == 1:
		return list(poly)
	r0 = poly_div(poly, tree[len(tree)-2-depth][2*index])
	r1 = poly_div(poly, tree[len(tree)-2-depth][2*index+1])
	list1 = mult_eval_recurse(r0, points[:num//2], tree, depth + 1, 2*index)
	list2 = mult_eval_recurse(r1, points[num//2:], tree, depth + 1, 2*index+1)
	return (list1 + list2)[:num] 
   

