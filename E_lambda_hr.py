from harmonic_oscillator import *
import matplotlib.pyplot as plt

data = np.arange(-6, 6, 0.1)

N = 20

def sq(x):
	return x**2


pot = sq(data)
plt.plot(data, pot)

H = matrix(N, 0)
M, vectors, counter = QR(H, 10**(-6))
eig_v_0 = [0 for _ in range(N)]
for i in range(N):
	eig_v_0[i] = M[i, i]*5
eig_v_0.sort()


def plot_hr(n, lamb):
	shift = 0
	H = matrix(N, lamb)
	M, vectors, counter = QR(H, 10**(-6))
	eig_v_computed = [0 for _ in range(N)]
	for i in range(N):
		eig_v_computed[i] = M[i, i]
	eig_v_computed.sort()
	for i in range(n):
		if i == 0:
			y = np.sqrt(eig_v_computed[i]*5)
			plt.hlines(eig_v_computed[i]*5, -1*y, y, color = col,label=str(lamb))
		else:
			y = np.sqrt(eig_v_computed[i]*5)
			plt.hlines(eig_v_computed[i]*5, -1*y, y, color = col)
		if lamb != 0 and i > 0:
			if col == 'g':
				shift = 0.3
			plt.arrow((-1.5+i/2+shift)/2, eig_v_0[i], 0, eig_v_computed[i]*5-eig_v_0[i], color = col, length_includes_head=True, head_length=.5, head_width=.2)

plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,
    labelbottom=False)
col = 'k'
plot_hr(6, 0)
col = 'r'
plot_hr(6, 0.02)
col = 'g'
plot_hr(6, 0.1)
plt.legend(title=r'$\lambda$')
plt.ylim(0, 30)
plt.yticks(eig_v_0[0:6], [r'$E_%i$' %i for i in range(6)])
plt.xlim(-6, 6)
plt.show()

