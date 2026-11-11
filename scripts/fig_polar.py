# C. Inzunza, pseudostress-based fully mixed FEM for stationary Biot-PNP (submitted), Example 3:
# ionic polarization in a charged porous medium, V = 0.9. Recomputed here in primal form (P1 Newton for
# Poisson-Nernst-Planck with the Darcy advection neglected, then P2-P1 Biot driven by the electric terms).
import numpy as np, sys, triangle, scipy.sparse as sp
sys.path.insert(0, 'scripts'); from sitefig import *
from skfem import *
from skfem.helpers import dot, grad, div, sym_grad, ddot
mu, lam, alpha, c0, delta, kap, nu_, k1, k2, eps = 1, 2, .1, .1, .1, 1, 1, 1, .7, 2e-4
Cf, R, eta, fsrc, xi_b = 0.04, 0.2, 0.02, 0.02, 0.02
T = triangle.triangulate({'vertices': np.array([(0, 0), (1, 0), (1, 1), (0, 1)]), 'segments': [[0, 1], [1, 2], [2, 3], [3, 0]]}, 'pq30a0.0002')
m = MeshTri(T['vertices'].T, T['triangles'].T); print('triangles', m.t.shape[1])
b = Basis(m, ElementTriP1(), intorder=4)
fchi = lambda x, y: -Cf/2*(1 + np.tanh((R**2 - (x - .5)**2 - (y - .5)**2)/(2*R*eta)))
Mm = asm(BilinearForm(lambda u, v, w: u*v), b); Kk = asm(BilinearForm(lambda u, v, w: dot(grad(u), grad(v))), b)
Fc = asm(LinearForm(lambda v, w: fchi(w.x[0], w.x[1])*v), b)
Fs = asm(LinearForm(lambda v, w: fsrc*v), b)
@BilinearForm
def drift(u, v, w): return u*dot(w['gchi'], grad(v))          # d/dxi of xi grad(chi).grad(w)
@BilinearForm
def dchi(u, v, w): return w['xi']*dot(grad(u), grad(v))         # d/dchi
N = b.N; bd = b.get_dofs().all()
def solve_pnp(V, Cscale, x0):
    chi, x1, x2 = x0
    for it in range(40):
        gchi = b.interpolate(chi).grad
        Dr = asm(drift, b, gchi=gchi)
        Rchi = eps*Kk @ chi - Mm @ (x1 - x2) - Cscale*Fc
        R1 = Mm @ x1 + k1*(Kk @ x1 + Dr @ x1) - Fs
        R2 = Mm @ x2 + k2*(Kk @ x2 - Dr @ x2) - Fs
        J = sp.bmat([[eps*Kk, -Mm, Mm],
                     [k1*asm(dchi, b, xi=b.interpolate(x1)), Mm + k1*(Kk + Dr), None],
                     [-k2*asm(dchi, b, xi=b.interpolate(x2)), None, Mm + k2*(Kk - Dr)]], 'csr')
        Rv = np.r_[Rchi, R1, R2]
        D = np.r_[bd, bd + N, bd + 2*N]
        dx = solve(*condense(J, -Rv, x=np.zeros(3*N), D=D))
        chi, x1, x2 = chi + dx[:N], x1 + dx[N:2*N], x2 + dx[2*N:]
        if np.linalg.norm(dx)/np.linalg.norm(np.r_[chi, x1, x2]) < 1e-10: break
    return (chi, x1, x2), it
X = b.doflocs
V = 0.9
state = (V*(X[0] - .5)*0, xi_b*np.ones(N), xi_b*np.ones(N))
for s in np.linspace(1/6, 1, 6):            # load continuation on the fixed charge and the voltage
    st = list(state); st[0] = st[0].copy(); st[0][bd] = s*V*(X[0, bd] - .5)
    state, it = solve_pnp(s*V, s, st); print('step', round(s, 3), 'newton', it)
chi, x1, x2 = state
print('chi range', chi.min(), chi.max(), 'charge diff', (x1 - x2).min(), (x1 - x2).max())
# ---- mechanics: steady Biot driven by Maxwell and osmotic stresses ----
ev, ep = ElementVector(ElementTriP2()), ElementTriP1()
bu, bp = Basis(m, ev, intorder=4), Basis(m, ep, intorder=4)
K = asm(BilinearForm(lambda u, v, w: 2*mu*ddot(sym_grad(u), sym_grad(v)) + lam*div(u)*div(v)), bu)
C = asm(BilinearForm(lambda u, q, w: alpha*div(u)*q), bu, bp)
Ap = asm(BilinearForm(lambda p, q, w: kap/nu_*dot(grad(p), grad(q)) + c0*p*q), bp)
gc = bu.interpolate(bu.project(lambda x: 0*x) if False else np.zeros(bu.N))  # placeholder
chi_u = b.interpolate(chi); dq = b.interpolate(x1 - x2)
# evaluate fields at the vector basis quadrature (same mesh, same intorder)
@LinearForm
def Fel(v, w):
    g = w['g']; s = w['s']
    tau = eps*np.einsum('i...,j...->ij...', g, g); trp = 0.5*eps*(g[0]**2 + g[1]**2) + delta*s
    return -(ddot(tau, grad(v)) - trp*div(v))
F = asm(Fel, bu, g=chi_u.grad, s=dq.value)
A = sp.bmat([[K, -C.T], [C, Ap]], 'csr')
x = solve(*condense(A, np.r_[F, np.zeros(bp.N)], D=bu.get_dofs().all()))
u, p = x[:bu.N], x[bu.N:]
c = m.p[:, m.t].mean(1)
U = bu.probes(c) @ u; n = c.shape[1]; umag = np.hypot(U[:n], U[n:])
gp = bp.interpolate(p).grad  # Darcy flux z = -kappa/nu grad p, element mean
zmag = np.hypot(gp[0].mean(1), gp[1].mean(1))*kap/nu_
print('|u| max', umag.max(), '|z| max', zmag.max())
# plot on a coarser triangulation (fields sampled from the fine solution)
Tp = triangle.triangulate({'vertices': np.array([(0, 0), (1, 0), (1, 1), (0, 1)]), 'segments': [[0, 1], [1, 2], [2, 3], [3, 0]]}, 'pq30a0.0011')
mp = MeshTri(Tp['vertices'].T, Tp['triangles'].T); cp = mp.p[:, mp.t].mean(1); print('plot triangles', mp.t.shape[1])
el = m.element_finder()(cp[0], cp[1])
chic = b.probes(cp) @ chi
U = bu.probes(cp) @ u; n = cp.shape[1]; umag = np.hypot(U[:n], U[n:]); zmag = zmag[el]
m = mp
W, gap = 190, 24; body = ''
for k, (v, pal, lab, kw) in enumerate([(chic, 'b', 'χ', dict()), (umag, 's', '|u|', dict(lo=0)), (zmag, 's', '|z|', dict(lo=0))]):
    P = Panel((0, 1, 0, 1), W, k*(W + gap), 0)
    body += tri_paths(P, m.p, m.t, bands(v, **kw), pal)
    body += label(k*(W + gap) + W/2, W + 22, lab)
r = R*W
for k in range(3):
    body += '<circle class="contour" cx="%s" cy="%s" r="%s"/>' % (fmt(k*(W + gap) + W/2), fmt(W/2), fmt(r))
standalone('assets/img/figs/bpnp-polarization.svg', 3*W + 2*gap, W + 30, body)
