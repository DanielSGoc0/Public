# A program that alternates between taking self-convolutions and squaring a distribution on Z/ZN.

from tkinter import *
import tkinter
import copy
import random
import math
import numpy
import scipy

# N = 16
# N = 32
N = 2048 + 1
mu = [0.0 for _ in range(N)]
# mu = [0.17054627, 0.01737024, 0.00637775, 0.15565203, 0.03952459, 0.002047, 0.11832952, 0.07493073, 0.00099]
var = 10000

def normal_mu():
	global N, mu, var

	for i in range(N):
		mu[i] = math.exp(-i*i/var)

def poisson_mu():
	global N, mu, var
	r = 0.9

	for i in range(N):
		theta = (math.pi * i)/(N - 1)
		mu[i] = (1 - r*r)/(1 - 2*r*math.cos(theta) + r*r)


def normal_mu2():
	global N, mu, var

	for i in range(N):
		mu[i] = math.exp(-i*i/var) - 2*math.exp(-i*i/var * 1.1) + 4*math.exp(-i*i/var * 1.3)
		# mu[i] = math.exp(-i*i/var)

def laplace_mu():
	global N, mu, var

	for i in range(N):
		mu[i] = math.exp(-i/N * 300)

def other_mu():
	global N, mu

	mu_transform = [0.0 for i in range(N)]
	c = 1

	for i in range(N):
		# mu_transform[i] = 1.0/(1000.0 + i*i) * math.exp(-i*i/1000)
		mu_transform[i] = 100.0/math.cosh(i/10.0) + 1.0/math.cosh(i/20.0)

	mu = scipy.fft.idct(mu_transform, type=1)

def random_mu():
	global N, mu

	while True:
		mu = [random.random() for _ in range(N)]
		mu[N//4] = 0
		mu_transform = scipy.fft.dct(mu, type=1)
		minim = 1.0
		for j in range(N):
			minim = min(minim, mu_transform[j])

		if minim > 0.0:
			break	
	# mu[0] = mu[0] - minim + random.random()*3

def partial_convolution(a):
	global N, mu

	mu_transform = scipy.fft.dct(mu, type=1)
	# print(mu)
	# print(mu_transform)

	minim = 1.0
	for i in range(N):
		minim = min(minim, mu_transform[i])
		mu_transform[i] = max(0, mu_transform[i])**(1 + a)

	print("MINIM = ", minim)
	mu = scipy.fft.idct(mu_transform, type=1)	

def MGF(a):
	global N, mu
	res = 0.0

	for i in range(N):
		s = 1.0
		if i == 0 or i == N - 1:
			s = 1.0/2
		res = res + s * (math.exp(a * i) + math.exp(-a * i)) * mu[i]

	return res

def random_mu2():
	global N, mu
	mu = [random.random() for _ in range(N)]

def random_mu3():
	global N, mu, var

	for i in range(N):
		mu[i] = math.exp(-i*i/var) * (1.5 - random.random())

def step_h(h):
	global N, mu
	for i in range(N):
		mu[i] = max(0, mu[i])**(1 + h)

	mu_transform = scipy.fft.dct(mu, type=1)
	# print(mu)
	# print(mu_transform)

	for i in range(N):
		mu_transform[i] = max(0, mu_transform[i])**(1 + h)

	mu = scipy.fft.idct(mu_transform, type=1)

def normalize():
	global mu, N
	norm = 0.0
	for i in range(N):
		norm = norm + mu[i]
	norm = 2*norm
	norm = norm - mu[0] - mu[N - 1]
	# norm = numpy.linalg.norm(mu)
	for i in range(N):
		mu[i] = mu[i] / norm

def entropy():
	global mu, N
	res = 0
	for i in range(N):
		res = res - mu[i] * math.log(mu[i])
	return res

def run_steps(K, h):
	global mu, N

	# poisson_mu()
	normal_mu()
	normalize()
	print(mu[0], mu[N//4], mu[N//2], mu[N//4 * 3], mu[N - 1])
	mu_transform = scipy.fft.dct(mu, type=1)
	minim = 1.0
	for j in range(N):
		minim = min(minim, mu_transform[j])
	print(minim)
	# print(mu)
	# saved = []
	# for i in range(N):
	# 	saved.append(mu[i])
	# for i in range(10):
	# 	mu = []
	# 	for j in range(N):
	# 		mu.append(saved[j])

	# 	l =  0.1**i

	# 	mu_transform = scipy.fft.dct(mu, type=1)
	# 	for i in range(N):
	# 		mu_transform[i] = mu_transform[i]**(1 + l)
	# 	mu = scipy.fft.idct(mu_transform, type=1)
	# 	print(l)
	# 	print((mu[N//4] - saved[N//4])/l)



	# for i in range(K):
	# 	normalize()
	# 	print(mu)
	# 	# print(entropy())
	# 	step_h(h)

	# print(mu)

def print_out():
	global N, mu, var

	tup = []
	for i in range(10):
		tup.append(mu[i] * math.exp(i*i/var))
	print(tup)


def run_steps2(K):
	global mu, N

	random_mu2()
	# normal_mu()
	# random_mu3()
	# normal_mu2()
	# normalize()
	# print(mu)

	for i in range(K):
		print(mu)
		# print_out()
		# print(entropy())
		step_h(1)
		normalize()

	print(mu)

# run_steps2(200)


# mu_transform = scipy.fft.dct(mu, type=1)
# minim = 1.0
# print(mu_transform)
# for j in range(N):
# 	minim = min(minim, mu_transform[j])
# print(minim)

# run_steps(5, 1.0/1000)

# normal_mu2()
other_mu()
normalize()
old_MGF = [MGF(0.01 / 1.3**i) for i in range(0, 10)]
print(old_MGF)
a = -0.5
partial_convolution(a)

new_MGF = [MGF(0.01 / 1.3**i) for i in range(0, 10)]
print(new_MGF)
for i in range(0, 10):
	print(old_MGF[i], old_MGF[i]**(1 + a) - new_MGF[i])

# normal_mu()
# print(mu)
# mu_transform = scipy.fft.dct(mu, type=1)
# minim = 1.0
# for j in range(N):
# 	minim = min(minim, mu_transform[j])
# print(minim)

# print(scipy.fft.dct(mu, type=1))
