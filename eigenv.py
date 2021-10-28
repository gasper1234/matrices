from harmonic_oscillator import *
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm
#from h_map import *

N = 200
H = matrix(N, 0.2)
a, b = np.linalg.eig(H)
#M, vectors, counter = QR(H, 10**(-6))
M = np.array([[0. for _ in range(N)] for _ in range(N)])
for i in range(N):
	M[i][i] = a[i]
ord = order(M)
A = normalize(b)

h_map = np.array([[0. for _ in range(N)] for _ in range(N)])

for i in range(1, N+1):
	for j in range(N):
		h_map[N-i, j] = A[j][ord[i-1]]


h_map = h_map[140:]

for i in range(len(h_map)):
	for j in range(len(h_map[0])):
		h_map[i][j] = np.sqrt(h_map[i][j])


fig, ax = plt.subplots()
tick_loc = [i-1 for i in range(len(h_map), 0, -5)]
plt.yticks(tick_loc, [i*5 for i in range(len(tick_loc))])
cax = ax.imshow(h_map, cmap='inferno')
ax.set_xlabel('l. stanja nemotenega HO')
ax.set_ylabel('l. stanja HO z motnjo')
cbar = fig.colorbar(cax, orientation="horizontal", ticks = [round(np.sqrt(i/10), 2) for i in range(0, 11, 2)])
cbar.set_ticklabels([str(i/10) for i in range(0, 11, 2)])
ax.set_title(r'$c_n^2$', size=10)
plt.show()
'''
fig, ax = plt.subplot()
tick_loc = [i-1 for i in range(len(h_map), 0, -5)]
plt.yticks(tick_loc, [i*5 for i in range(len(tick_loc))])
cax = ax.imshow(h_map, cmap='inferno')
plt.xlabel('l. stanja nemotenega HO')
plt.ylabel('l. stanja HO z motnjo')
divider = make_axes_locatable(ax)
divider.append_axes("top", size="10%", pad=0.5)
fig.colorbar(cax, orientation="horizontal", ticks = [round(np.sqrt(i/10), 2) for i in range(0, 11, 2)])
plt.title(r'$c_n^2$', size=10)
plt.show()
'''