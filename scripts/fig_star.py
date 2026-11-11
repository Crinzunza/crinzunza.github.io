# Biot-PNP (C. Inzunza, submitted), Example 2: internal layer on a three-quarter star.
# Exact fields from the manufactured solution; mesh refined towards the layer.
import numpy as np, triangle, sys
sys.path.insert(0, 'scripts'); from sitefig import *
from skfem import MeshTri
a, A, w = 1e-3, 5e-3, 0.1
V = [(0.0, 0.0)] + [((1.45 if j % 2 == 0 else 0.95)*np.cos(j*np.pi/8), (1.45 if j % 2 == 0 else 0.95)*np.sin(j*np.pi/8)) for j in range(13)]
V = np.array(V); seg = np.c_[np.arange(len(V)), np.roll(np.arange(len(V)), -1)]
T = triangle.triangulate({'vertices': V, 'segments': seg}, 'pq30a0.08')
m = MeshTri(T['vertices'].T, T['triangles'].T)
H = lambda x, y: np.tanh((x + y)/w)
def charge(c): return 2*A*H(*c)
def dispmag(c):   # |phi| = eps |grad chi|, chi = a(1 - H^2), eps = 1
    h = H(*c); return a*2*np.abs(h)*(1 - h**2)*np.sqrt(2)/w
for it in range(9):
    c = m.p[:, m.t].mean(1)
    e1 = m.p[:, m.t[1]] - m.p[:, m.t[0]]; e2 = m.p[:, m.t[2]] - m.p[:, m.t[0]]
    hT = np.sqrt(np.abs(e1[0]*e2[1] - e1[1]*e2[0]))
    ind = hT**2*(1 - H(*c)**2)/w            # |grad H| scaled by element area
    mark = np.nonzero(ind >= 0.8*ind.mean())[0]
    if m.t.shape[1] > 1100: break
    m = m.refined(mark)
print('triangles', m.t.shape[1])
c = m.p[:, m.t].mean(1)
W = 330; view = (-1.5, 1.5, -1.5, 1.5); gap = 40
P1 = Panel(view, W, 0, 0); P2 = Panel(view, W, W + gap, 0)
body = tri_paths(P1, m.p, m.t, bands(charge(c), sym=True), 'b')
body += tri_paths(P2, m.p, m.t, bands(dispmag(c), lo=0), 's')
body += label(W/2, W + 22, 'ξ₁ − ξ₂') + label(W + gap + W/2, W + 22, '|φ|')
standalone('assets/img/figs/bpnp-star.svg', 2*W + gap, W + 30, body)
