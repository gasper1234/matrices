from harmonic_oscillator import *
import matplotlib.pyplot as plt


for N in range(10, 100, 10):
	H = matrix(N, 0.1)
	#eig_v_np, b = np.linalg.eig(H)
	M, vectors, counter = QR(H, 10**(-6))
	eig_v_computed = [0 for _ in range(N)]


	for i in range(N):
		eig_v_computed[i] = M[i, i]

	eig_v_computed.sort()

	N_l = [i for i in range(N)]

	plt.plot(N_l, eig_v_computed, '-', label=str(N))

plt.xlabel('N')
plt.ylabel(r'$E_n$')
plt.legend()
plt.ylim(0, 120)
plt.xlim(0, 40)
plt.legend(title=r'dim')
plt.show()