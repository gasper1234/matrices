import numpy as np
from math import factorial as fac


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
	return 1 / 2**4 * np.sqrt((2**i * fac(i)) / (2**j * fac(j))) * (a + b + c + d + e)

def matrix(N, lamb):
	A = np.array([[0. for _ in range(N)] for _ in range(N)])
	for i in range(N):
		for j in range(N):
			sum = H_0(i, j) + lamb * q_4(i, j)
			A[i][j] = sum
	return A


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

def sum_nondiag(M):
	N = len(M)
	sum = 0
	for i in range(N):
		for j in range(N):
			if i != j:
				sum += M[i][j]
	return sum

def QR(M):
	A = np.copy(M)
	N = len(A)
	for i in range(110):
		Q = np.identity(N)
		for i in range(len(A)-1):
			Q = np.dot(Housh_matrix(A, i), Q)
		A = np.dot(Q, np.dot(A, Q.transpose()))
#		print(sum_nondiag(A))
	return A

def diag(M):
	for i in range(len(M)):
		print(M[i][i])

H = matrix(150, 0.1)

A = H - 1.53564828 * np.identity(len(H))

print(np.linalg.det(A))

a, b = np.linalg.eig(H)

print(a)

diag(QR(H))
