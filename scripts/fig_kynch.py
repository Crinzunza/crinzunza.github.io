# Careaga & Inzunza, heat-assisted sedimentation (arXiv:2610.28903), Example 2, batch Kynch test:
# vessel (a), phi0 = 0.3, fbk = v0 phi (1-phi)^2, D1 = 0.015 exp(-phi^2), dt = 0.005.
# Recomputed here for the solids equation alone (no bulk flow, no heat coupling), P1 elements.
import numpy as np, sys, triangle
sys.path.insert(0, 'scripts'); from sitefig import *
from skfem import *
from skfem.helpers import dot, grad
v0, dt = 0.5, 0.005
V = np.array([(-0.5, 0), (0.5, 0), (1, 2), (-1, 2)])
T = triangle.triangulate({'vertices': V, 'segments': [[0, 1], [1, 2], [2, 3], [3, 0]]}, 'pq30a0.003')
m = MeshTri(T['vertices'].T, T['triangles'].T); print('triangles', m.t.shape[1])
b = Basis(m, ElementTriP1())
M = asm(BilinearForm(lambda u, v, w: u*v), b)
phi = 0.3*np.ones(b.N)
@BilinearForm
def Af(u, v, w): return w['D']*dot(grad(u), grad(v))
@LinearForm
def Ff(v, w):
    ph = w["ph"]; return -v0*ph*(1 - ph)**2*grad(v)[1]
times = [0.5, 1.5, 3.0, 5.0]; snaps = {}
for n in range(1, 1001):
    ph = b.interpolate(phi)
    A = asm(Af, b, D=0.015*np.exp(-ph**2))
    rhs = M @ phi + dt*asm(Ff, b, ph=ph)
    phi = solve(M + dt*A, rhs)
    t = round(n*dt, 6)
    if t in times: snaps[t] = phi.copy(); print(t, phi.min(), phi.max(), 'mass', (M @ phi).sum())
c = m.p[:, m.t].mean(1)
W, gap = 140, 18; body = ''
view = (-1.02, 1.02, -0.02, 2.02)
for k, t in enumerate(times):
    P = Panel(view, W, k*(W + gap), 0)
    vals = b.probes(c) @ snaps[t]
    body += tri_paths(P, m.p, m.t, bands(vals, lo=0, hi=0.9), 's')
    body += label(k*(W + gap) + W/2, P.h + 22, 't = %g s' % t)
standalone('assets/img/figs/kynch-batch.svg', 4*W + 3*gap, P.h + 30, body)
