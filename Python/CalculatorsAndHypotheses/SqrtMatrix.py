# Check if matrix:
# ((A, A#B, A#C), (B#A, B, B#C), (C#A, C#B, C))
# is positive semi-definite
import numpy as np
from scipy.stats import ortho_group
from scipy.linalg import sqrtm


def Gmean(A, B):
	Asqrt = sqrtm(A)
	Asqrt = (Asqrt + Asqrt.T) / 2.0
	Asqrtinv = np.linalg.inv(Asqrt)
	Asqrtinv = (Asqrtinv + Asqrtinv.T) / 2.0

	M = Asqrtinv @ B @ Asqrtinv
	M = sqrtm(M)
	M = (M + M.T) / 2.0
	M = Asqrt @ M @ Asqrt
	return M

N = 4
for attempt in range(10000):
	A = np.diag(np.random.random((N)))
	for i in range(N):
		A[i][i] = A[i][i]/(1.0 - A[i][i])
	U = ortho_group.rvs(dim=N)
	A = U @ A @ U.T

	B = np.diag(np.random.random((N)))
	for i in range(N):
		B[i][i] = B[i][i]/(1.0 - B[i][i])
	U = ortho_group.rvs(dim=N)
	B = U @ B @ U.T

	C = np.diag(np.random.random((N)))
	for i in range(N):
		C[i][i] = C[i][i]/(1.0 - C[i][i])
	U = ortho_group.rvs(dim=N)
	C = U @ C @ U.T

	M = np.zeros((3*N, 3*N), dtype='double')
	M[0:N, 0:N] = A
	M[0:N, N:(2*N)] = Gmean(A, B)
	M[0:N, (2*N):(3*N)] = Gmean(A, C)
	M[N:(2*N), 0:N] = M[0:N, N:(2*N)].T
	M[N:(2*N), N:(2*N)] = B
	M[N:(2*N), (2*N):(3*N)] = Gmean(B, C)
	M[(2*N):(3*N), 0:N] = M[0:N, (2*N):(3*N)].T
	M[(2*N):(3*N), N:(2*N)] = M[N:(2*N), (2*N):(3*N)].T
	M[(2*N):(3*N), (2*N):(3*N)] = C
	

	# evals = np.linalg.eigvalsh(M[0:(2*N), 0:(2*N)])
	evals = np.linalg.eigvalsh(M)
	print(evals)
	break


