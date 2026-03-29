# Checks if the two operations on (F, T) coordinate system are commutative.
# This coordinate system arises when analyzing the class of G(F, A) functions.

# It turns out that the two operations are commutative! (not that it matters much)

import numpy as np
from numpy import *


def getFA(F, T):
	return (F * np.exp(-T / 2.0) * (1 + T), 1 - F**2 * np.exp(-T) * T)

def getFT(F, A):
	v = np.sqrt(1.0 - 4*(1 - A)/F**2)
	T = (1.0 - v)/(1.0 + v)
	F = F * np.exp(T/2.0)/(1.0 + T)
	return (F, T)

def transformFT(F, T, t):
	return (F, T + t)

def transformFA(F, A, s):
	return (np.exp(-s / 2.0) * F, 1.0 - np.exp(-s) * (1 - A))


F = 0.8
T = 0.9
t = 0.8
s = 0.9

(F0, A0) = getFA(F, T)
(F1, A1) = getFA(F, T + t)
(F2, A2) = transformFA(F1, A1, s)
(F1p, A1p) = transformFA(F0, A0, s)
(F1p, T1p) = getFT(F1p, A1p)
T1p += t
(F2p, A2p) = getFA(F1p, T1p)

print(F2, F2p)
print(A2, A2p)
