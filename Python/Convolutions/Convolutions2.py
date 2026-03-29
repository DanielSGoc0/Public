# A program that alternates between taking self-convolutions and squaring a distribution on Z/ZN.


from tkinter import *
import tkinter
import copy
import random
import math
import numpy as np
import scipy

N = 2**10
tab = [0.0 for _ in range(N)]
gamma = -1.0
sigma2 = 0.0001
h = 0.1

sigma_true = sigma2 / h**2
power_true = (1 + sigma_true)**(-gamma)

weights = [math.exp(-min(i, N - i)**2 * sigma_true / 2) for i in range(N)]
tab = [random.lognormvariate(0, 1) for i in range(N)]


# print(weights)


def convolution():
	global N, tab, weights

	tab = scipy.fft.fft(tab)
	for i in range(N):
		tab[i] *= weights[i]
	tab = scipy.fft.ifft(tab)
	tab = tab.real

def raise_to_power():
	global N, tab, power_true

	for i in range(N):
		tab[i] = max(0, tab[i])**power_true

def normalize():
	global tab, N
	norm = np.linalg.norm(tab)
	for i in range(N):
		tab[i] = tab[i] / norm

def run_steps(K):
	global tab
	for i in range(K):
		raise_to_power()
		convolution()
		normalize()

		print(tab[0:100])

# print(tab[0:100])
run_steps(1000000)
