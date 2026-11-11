# Convective Brinkman-Forchheimer flow in the non-convex channel of
# Caucao, Gatica, Gatica & Inzunza (CI2MA pp26-12, Example 4), recomputed with
# Taylor-Hood P2-P1 elements, Picard iteration and gradient-jump adaptivity.
import numpy as np, pickle
from skfem import *
from skfem.helpers import dot, grad, div, ddot
import scipy.sparse as sp

mu = 0.01
def poro(x, y):
    re = 0.65 - np.exp(-50*(x-1.2)**2 - 50*(y-0.25)**2) - 0.9*np.exp(-70*(x-1.85)**2 - 20*(y-0.10)**2)
    return re + (1 - re)*np.exp(y - 1)
Dc = lambda r: 150*((1-r)/r)**2
Fc = lambda r: 1.75*((1-r)/r)

def base_mesh():
    m = MeshTri.init_tensor(np.linspace(0, 2.5, 51), np.linspace(0, 0.3, 7))
    c = m.p[:, m.t].mean(1)
    cav = (c[0] > 0.6) & (c[0] < 0.8) & (c[1] < 0.2)
    return m.remove_elements(np.nonzero(cav)[0])

def solve(m, w_prev=None):
    eu = ElementVector(ElementTriP2()); ep = ElementTriP1()
    bu = Basis(m, eu, intorder=4); bp = Basis(m, ep, intorder=4)
    xu = bu.doflocs
    @BilinearForm
    def A(u, v, w):
        r = poro(w.x[0], w.x[1])
        wv = w['wv']; nw = np.sqrt(wv[0]**2 + wv[1]**2)
        conv = np.einsum('i...,ji...->j...', wv, grad(u))
        return mu*ddot(grad(u), grad(v)) + dot(conv, v) + (Dc(r) + Fc(r)*nw)*dot(u, v)
    @BilinearForm
    def B(u, q, w):
        return -div(u)*q
    # boundary dofs
    D = bu.get_dofs().all()
    def inflow(x, y): return 10*y*(0.3-y)
    uD = np.zeros(bu.N)
    lr = bu.get_dofs(lambda x: (x[0] < 1e-9) | (x[0] > 2.5-1e-9))
    ux = lr.nodal['u^1'] if 'u^1' in lr.nodal else None
    for key in ('u^1',):
        pass
    xs = xu[0]; ys = xu[1]
    # x-components are even dofs in ElementVector ordering? use nodal_dofs via get_dofs
    dofs_lr = bu.get_dofs(lambda x: (x[0] < 1e-9) | (x[0] > 2.5-1e-9))
    ix = np.concatenate([dofs_lr.nodal.get('u^1', np.array([], int)), dofs_lr.facet.get('u^1', np.array([], int))])
    uD[ix] = inflow(xs[ix], ys[ix])
    w = np.zeros(bu.N) if w_prev is None else w_prev
    Bm = asm(B, bu, bp)
    for it in range(30):
        Am = asm(A, bu, wv=bu.interpolate(w))
        K = sp.bmat([[Am, Bm.T], [Bm, 1e-10*asm(BilinearForm(lambda p, q, _: p*q), bp)]], 'csr')
        rhs = np.zeros(K.shape[0]); x0 = np.zeros(K.shape[0]); x0[:bu.N] = uD
        sol = solve_linear = None
        Dall = np.concatenate([D, [bu.N]])  # pin one pressure dof
        sol = solve_(K, rhs, x0, Dall)
        unew = sol[:bu.N]
        err = np.linalg.norm(unew - w)/max(np.linalg.norm(unew), 1e-14)
        w = unew
        if err < 1e-7: break
    return bu, bp, w, sol[bu.N:], it

def solve_(K, rhs, x0, D):
    from skfem import condense, solve as fsolve
    return fsolve(*condense(K, rhs, x=x0, D=D))

def estimator(m, bu, u):
    # gradient-jump indicator on interior facets
    fb = [InteriorFacetBasis(m, bu.elem, side=i, intorder=4) for i in (0, 1)]
    @Functional
    def jump(w):
        h = w.h; n = w.n
        g = grad(w['u0']) - grad(w['u1'])
        j = np.einsum('ij...,j...->i...', g, n)
        return h*dot(j, j)
    eta_f = jump.elemental(fb[0], u0=fb[0].interpolate(u), u1=fb[1].interpolate(u))
    eta = np.zeros(m.t.shape[1])
    f2t = m.f2t
    facets = fb[0].find
    for k, f in enumerate(facets):
        for t in f2t[:, f]:
            if t >= 0: eta[t] += 0.5*eta_f[k]
    # element residual of the drag terms (scaled), to pick up the low-porosity patches
    @Functional
    def res(w):
        r = poro(w.x[0], w.x[1]); uu = w['u']; nu = np.sqrt(uu[0]**2+uu[1]**2)
        return w.h**2*(Dc(r)+Fc(r)*nu)**2*dot(uu, uu)*1e-2
    eta += res.elemental(Basis(m, bu.elem, intorder=4), u=bu.interpolate(u))
    return np.sqrt(eta)

def dorfler(eta, theta=0.5):
    idx = np.argsort(-eta); c = np.cumsum(eta[idx]**2)
    return idx[:np.searchsorted(c, theta*c[-1]) + 1]

m = base_mesh(); w = None; levels = []
for lvl in range(9):
    bu, bp, u, p, it = solve(m)
    eta = estimator(m, bu, u)
    levels.append(dict(p=m.p.copy(), t=m.t.copy(), u=bu.interpolate(u), uvals=u, pres=p))
    vel = np.sqrt(sum(c**2 for c in bu.probes(m.p[:, m.t].mean(1)) @ u for c in [] ) ) if False else None
    print(lvl, m.t.shape[1], 'picard', it, 'eta', np.sqrt((eta**2).sum()))
    if m.t.shape[1] > 9000: break
    m = m.refined(dorfler(eta, 0.4))
# element-centroid velocity magnitude and pressure for each level
out = []
for L in levels:
    mm = MeshTri(L['p'], L['t'])
    bu = Basis(mm, ElementVector(ElementTriP2())); bp = Basis(mm, ElementTriP1())
    c = mm.p[:, mm.t].mean(1)
    U = bu.probes(c) @ L['uvals']; n = c.shape[1]
    vm = np.sqrt(U[:n]**2 + U[n:]**2) if U.shape[0] == 2*n else None
    out.append(dict(p=L['p'], t=L['t'], vel=vm, pres=bp.probes(c) @ L['pres'], rho=poro(*c)))
pickle.dump(out, open('scripts/channel.pkl', 'wb'))
