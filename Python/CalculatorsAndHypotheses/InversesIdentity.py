# Used to play around with certain algebraic identity
import numpy as np
import quaternion
import random
# from numpy.linalg import inv

N = 4

# A = np.array([[random.random() for j in range(N)] for i in range(N)])
A = np.array([random.random() for i in range(N)])
# A = quaternion.quaternion(random.random(), random.random(), random.random(), random.random())
# A = A + np.transpose(A)
# B = np.array([[random.random() for j in range(N)] for i in range(N)])
B = np.array([random.random() for i in range(N)])
# B = quaternion.quaternion(random.random(), random.random(), random.random(), random.random())
# B = B + np.transpose(B)
# B = np.eye(N)
# B = 2*B

def inv(a):
	len = a[0]**2 + a[1]**2 + a[2]**2 + a[3]**2
	# return np.array([a[0]/len, a[1]/len, a[2]/len, a[3]/len])
	return np.array([a[0]/len, -a[1]/len, -a[2]/len, -a[3]/len])


X = inv(A + B) + inv(inv(A) + inv(B))
Y = inv(inv(A + inv(B)) + inv(inv(A) + B))
print(X)
print(Y)
print(X[0] / Y[0], X[1] / Y[1], X[2] / Y[2], X[3] / Y[3])
# print()
# print(np.dot(X, Y))
# print(np.dot(Y, X))
# print(np.dot(X, Y) - np.dot(Y, X))
