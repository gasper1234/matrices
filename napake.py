from harmonic_oscillator import *
import matplotlib.pyplot as plt


N = 80
H = matrix(N, 0.1)
eig_v_np, b = np.linalg.eig(H)

reser = np.array([[0. for _ in range(N)] for _ in range(N)])
for i in range(N):
	reser[i][i] = eig_v_np[i]
ord = order(reser)

reser = np.copy(eig_v_np)
for i in range(N):
	eig_v_np[i] = reser[ord[i]]

for k in range(2, 10, 4):
	M, vectors, counter = QR(H, 10**(-k))
	eig_v_computed = [0 for _ in range(N)]

	ord = order(M)

	cop = []
	for i in range(N):
		cop.append(M[i][i])
	reser = np.copy(eig_v_np)
	for i in range(N):
		eig_v_computed[i] = cop[ord[i]]

	#error
	error_abs = abs(eig_v_np-eig_v_computed)
	error_rel = abs(eig_v_np-eig_v_computed) / eig_v_np
	N_l = [i for i in range(N)]

	plt.plot(N_l, error_rel, '-', label='rel. natančnost za '+r'$10^{-%s}$' % (k))
	plt.plot(N_l, error_abs, '-', label='abs. natančnost za '+r'$10^{-%s}$' % (k))

plt.xlabel('N-ta l. vrednost')
plt.ylabel('napaka')
plt.legend()
plt.yscale('log')
plt.ylim(10 ** (-17), 10 ** (-2))
plt.show()

'''
	for i in range(N):
		eig_v_computed[i] = M[i, i]

	eig_v_compare = np.copy(eig_v_computed)

	for i in range(N):
		for j in range(N):
			if abs(eig_v_np[i] - eig_v_computed[j]) < abs(eig_v_np[i] - eig_v_compare[i]):
				eig_v_compare[i] = eig_v_computed[j]

'''