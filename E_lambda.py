from harmonic_oscillator import *
import matplotlib.pyplot as plt

N = 50

H = matrix(N, 0)
eig_v_0 = [0 for _ in range(N)]
for i in range(N):
	eig_v_0[i] = H[i, i]
eig_v_0.sort()

print(eig_v_0)

lambs = [0, 0.01, 0.2]
colors = ['blue', 'green', 'red']

for l in range(len(lambs)):
	lamb = lambs[l]
	for N in range(30, 51, 10):
		if lamb == 0 and N < 50:
			continue
		H = matrix(N, lamb)
		#eig_v_np, b = np.linalg.eig(H)
		M, vectors, counter = QR(H, 10**(-6))
		eig_v_computed = [0 for _ in range(N)]

		for i in range(N):
			eig_v_computed[i] = M[i, i]

		eig_v_computed.sort()

		#for i in range(N):
		#	eig_v_computed[i] -= eig_v_0[i]


		N_l = [i for i in range(N)]

		if N == 30:
			line = 'dotted'
		if N == 40:
			line = 'dashdot'
		if N == 50:
			line = 'dashed'

		plt.plot(N_l, eig_v_computed, color=colors[l], linestyle=line, label=str(lamb)+', '+str(N))

plt.xlabel('N')
plt.ylabel(r'$E_n$')
plt.legend(title=r'$\lambda, dim$')
plt.ylim(0, 100)
plt.show()