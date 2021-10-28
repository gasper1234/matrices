import numpy as np
from scipy.special import hermite
#from math import factorial as fac

def H_n(x, n):
	eks = np.e ** (-1/2 * x**2)
	if n != 0:
		return eks * hermite(n)(x) / n * (7-n) / 3
	else:
		return eks * hermite(n)(x) * (7-n) / 3

def fac(i, j):
	mult = 1
	if i > j:
		for k in range(i-j):
			mult *= i-k
	elif j > i:
		for k in range(j-i):
			mult *= j-k
		mult = 1/mult
	return mult


def H_0(i, j):
	a = 0
	if i == j:
		a = i + 1/2
	return a

def q_1(i, j):
	if abs(i-j) == 1:
		return 1/2 * np.sqrt(i+j+1)
	else:
		return 0

def q_2(i, j):
	a, b, c = 0, 0, 0
	if i == j:
		b = 2*j + 1
	if i == (j-2):
		a = np.sqrt(j * (j-1))
	if i == (j+2):
		c = np.sqrt((j+1) * (j+2))
	return 1/2 * (a+b+c)


def q_4(i, j):
	a, b, c, d, e = 0, 0, 0, 0, 0
	if i == (j+4):
		a = 1
	if i == (j+2):
		b = 4 * (2*j+3)
	if i == j:
		c = 12 * (2*j**2 + 2*j + 1)
	if i == (j-2):
		d = 16*j * (2*j**2 -  3*j + 1)
	if i == (j-4):
		e = 16*j * (j**3 - 6*j**2 + 11*j - 6)
	big_int = fac(i, j)
	if a+b+c+d+e == 0:
		return 0
	return 1 / 2**4 * np.sqrt(2 ** (i-j) * float(fac(i, j))) * (a + b + c + d + e)

def matrix_mono(N,lamb):
	A = np.array([[0. for _ in range(N)] for _ in range(N)])
	L = np.sqrt(np.sqrt(lamb))
	for i in range(N):
		for j in range(N):
			A[i][j] = L * q_1(i, j)
	return A

def matrix_sq(N,lamb):
	A = np.array([[0. for _ in range(N)] for _ in range(N)])
	L = np.sqrt(lamb)
	for i in range(N):
		for j in range(N):
			A[i][j] = L * q_2(i, j)
	return A

def matrix_H(N):
	A = np.array([[0. for _ in range(N)] for _ in range(N)])
	for i in range(N):
		A[i][i] = H_0(i, i)
	return A

#vrne Hamiltonian nasega problema
def matrix(N, lamb):
	A = np.array([[0. for _ in range(N)] for _ in range(N)])
	for i in range(N):
		for j in range(N):
			sum = H_0(i, j) + lamb * q_4(i, j)
			A[i][j] = sum
	return A

#vrne housholderjevo matriko, ki iznici n-ti stolpec M (od diagonale navzdol)
def Housh_matrix(M, n):
	N = len(M)
	c = np.array([M[i][n] for i in range(n, N)])
	A = np.array([[0. for _ in range(N)] for _ in range(N)])

	c_2 = 0
	for i in c:
		c_2 += i**2

	c[0] += np.sign(c[0]) * np.sqrt(c_2)

	c_2 = 0
	for i in c:
		c_2 += i**2

	for i in range(N):
		for j in range(N):
			if i == j:
				A[i][j] += 1
			if i >= n and j >= n:
				A[i][j] -= 2 * c[i-n]*c[j-n] / c_2
	return A

#sesteje izvendiagonalne clene
def sum_nondiag_sq(M):
	N = len(M)
	sum = 0
	for i in range(N):
		for j in range(N):
			if i != j:
				sum += (M[i][j])**2
	return sum

#vrne diagonalno obliko matrike (s QR razcepom)
def QR(M, eps):
	A = np.copy(M)
	N = len(A)
	count = 0
	Q_eig = np.identity(N)
	while  sum_nondiag_sq(A) > eps:
		count += 1
		Q = np.identity(N)
		for i in range(len(A)-1):
			Q = np.dot(Housh_matrix(A, i), Q)
		A = np.dot(Q, np.dot(A, Q.transpose()))
		Q_eig = np.dot(Q_eig, Q.transpose())
	return A, Q_eig, count

#sprinta diagonale
def diag(M):
	for i in range(len(M)):
		print(M[i][i])

#normalizira matriko
def normalize(M):
	A = np.copy(M)
	for i in range(len(M)):
		norm = 0
		for j in range(len(M)):
			norm += A[i][j]**2
		norm = np.sqrt(norm)
		for j in range(len(M)):
			A[i][j] /= norm
			A[i][j] *= A[i][j]		
	return A

def normalize_1(M):
	A = np.copy(M)
	for i in range(len(M)):
		norm = 0
		for j in range(len(M)):
			norm += A[i][j]**2
		norm = np.sqrt(norm)
		for j in range(len(M)):
			A[i][j] /= norm	
	return A

#vrne vrstni red lastnih vrednosti matrike po velikosti
def order(M):
	N = len(M)
	val = [M[i][i] for i in range(N)]
	max_ind = 0
	for i in range(N):
		if val[i] > val[max_ind]:
			max_ind = i
	ind = [max_ind for _ in range(N)]
	for i in range(N):
		for j in range(N):
			if i == 0:
				if val[j] < val[ind[i]]:
					ind[i] = j
			else:
				if val[j] < val[ind[i]] and val[j] > val[ind[i-1]]:
					ind[i] = j
	return ind

#testi

'''
N = 5
H = matrix(N, 0.01)
H_2 = np.dot(matrix_sq(N, 0.01),matrix_sq(N, 0.01))
print(H_2+matrix_H(N))
print(H)
'''


'''
M, vectors, c = QR(H_2, 10**(-10))
diag(M)
M, vectors, c = QR(H, 10**(-10))
diag(M)

#a, b = np.linalg.eig(H)

M, eig_v, counter = QR(H, 10**(-10))

diag(M)
print(eig_v)
'''
