from harmonic_oscillator import *
import matplotlib.pyplot as plt

data = np.arange(-6, 6, 0.1)

N = 5
lamb = 0

def sq(x, lamb):
	return (x**2+x**4*lamb) * 2

def H_n_dist(x, n, M):
	koef = np.copy(M)
	sum = 0
	for i in range(len(koef[:,n])):
		sum += koef[:,n][i] * H_n(x, i)
	return sum

fig, axs = plt.subplots(1, 3)
lambdas = [0, 0.02, 1]

for r in range(3):

	lamb = lambdas[r]

	#izračuna koeficiente funkcij motenega HO
	H = matrix(N, lamb)
	M, vectors, counter = QR(H, 10**(-6))
	ord = order(M)
	A = normalize_1(vectors)
	koef = np.array([[0. for _ in range(N)] for _ in range(N)])

	for i in range(N):
		for j in range(N):
			koef[j, i] = A[j][ord[i]]

	pot = sq(data, lamb)
	y_ticks = []

	for i in range(5):
		data_add = [i*7 for _ in range(len(data))]
		H = H_n_dist(data, i, koef)
		H = H + data_add
		axs[r].plot(data, H)
		y_ticks.append(data_add[0]+H_n_dist(-5, i, koef))

	axs[r].plot(data, pot, 'k')

	axs[r].tick_params(
	    axis='x',          # changes apply to the x-axis
	    which='both',      # both major and minor ticks are affected
	    bottom=False,      # ticks along the bottom edge are off
	    top=False,
	    labelbottom=False)
	axs[r].set_xlabel(r'$\lambda = $'+str(lamb))
	axs[r].set_ylim(-5, 35)
	axs[r].set_xlim(-5, 5)
	axs[r].set_yticks(y_ticks)
	axs[r].set_yticklabels([r'$\psi_%i$' %i for i in range(len(y_ticks))])
plt.show()
