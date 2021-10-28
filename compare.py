from harmonic_oscillator import *
import matplotlib.pyplot as plt


N = 500
lamb = 0.1
H = matrix(N, lamb)
SQ = matrix_sq(N, lamb)
MO = matrix_mono(N, lamb)
H_0 = matrix_H(N)

eig_v_4, b = np.linalg.eig(H)
eig_v_2, b = np.linalg.eig(np.dot(SQ, SQ) + H_0)
eig_v_1, b = np.linalg.eig(np.dot(np.dot(MO, MO), np.dot(SQ, SQ)) + H_0)

l = [eig_v_4, eig_v_2, eig_v_1]

for r in range(3):
	M = np.array([[0. for _ in range(N)] for _ in range(N)])
	for i in range(N):
		M[i][i] = l[r][i]
	ord = order(M)

	cop = np.copy(l[r])

	for i in range(N):
		l[r][i] = cop[ord[i]]

print(l)

#error
error_abs_2 = abs(l[0]-l[1])
error_rel_2 = abs(l[0]-l[1]) / l[0]
error_abs_1 = abs(l[0]-l[2])
error_rel_1 = abs(l[0]-l[2]) / l[0]
N_l = [i for i in range(N)]


plt.plot(N_l, error_rel_2, '-', label=r'rel. razlika $q^4$ in $q^2$')
plt.plot(N_l, error_abs_2, '-', label=r'abs. razlika $q^4$ in $q^2$')

plt.plot(N_l, error_rel_1, '-', label=r'rel. razlika $q^4$ in $q$')
plt.plot(N_l, error_abs_1, '-', label=r'rel. razlika $q^4$ in $q$')

plt.xlabel('N-ta l. vrednost')
plt.ylabel('napaka')
plt.legend(loc='upper left')
plt.yscale('log')
plt.xlim(-10, 160)
plt.show()