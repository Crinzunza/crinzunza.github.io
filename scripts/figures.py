# Builds the site figures from the numerical examples in Cristian Inzunza's CI2MA preprints.
import numpy as np, pickle
from skfem import MeshTri
OUT = '_includes/'
NB = 7

def svg_from(p, t, val, pal, view, width, lo=None, hi=None, extra='', cls='mesh', groups=None, slice_=False, quant=False):
    """p: 2xN points (model coords), t: 3xM, val: M values. view=(x0,x1,y0,y1)."""
    x0, x1, y0, y1 = view
    s = width/(x1-x0); H = (y1-y0)*s
    def tr(P):  # flip y
        return np.c_[(P[0]-x0)*s, (y1-P[1])*s]
    def paths(p, t, val):
        if quant:
            edges = np.quantile(val, np.linspace(0, 1, NB+1)[1:-1]); q = np.searchsorted(edges, val)
        else:
            lo_ = np.nanmin(val) if lo is None else lo; hi_ = np.nanmax(val) if hi is None else hi
            q = np.clip(((val-lo_)/(hi_-lo_)*NB).astype(int), 0, NB-1)
        P = tr(p); out = {}
        for k in range(t.shape[1]):
            tri = P[t[:, k]]
            if tri[:, 0].max() < -2 or tri[:, 0].min() > width+2 or tri[:, 1].max() < -2 or tri[:, 1].min() > H+2: continue
            out.setdefault(int(q[k]), []).append("M%.1f %.1fL%.1f %.1fL%.1f %.1fZ" % tuple(tri.ravel()))
        return "".join('<path class="%s%d" d="%s"/>' % (pal, b, "".join(v).replace('.0', '')) for b, v in sorted(out.items()))
    par = ' preserveAspectRatio="xMidYMid slice"' if slice_ else ''
    head = '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d"%s class="%s" aria-hidden="true" focusable="false">' % (width, round(H), par, cls)
    if groups:
        body = "".join('<g class="lvl lvl%d">%s</g>' % (i, paths(*g)) for i, g in enumerate(groups))
    else:
        body = paths(p, t, val)
    return head + body + extra + '</svg>', tr

L = pickle.load(open('scripts/channel.pkl', 'rb'))
vmax = max(l['vel'].max() for l in L)
# ---- hero: channel, three refinement levels, velocity magnitude ----
sel = [L[0], L[4], L[-1]]
svg, _ = svg_from(None, None, None, 's', (0, 2.5, 0, 0.3), 1500, lo=0, hi=vmax,
                  groups=[(l['p'], l['t'], l['vel']) for l in sel], cls='mesh hero-mesh', quant=True)
open(OUT+'hero-mesh.svg', 'w').write(svg); print('hero', len(svg), [l['t'].shape[1] for l in sel])
# ---- flow: velocity near the cavity, final mesh ----
f = L[-1]
svg, _ = svg_from(f['p'], f['t'], f['vel'], 's', (0.35, 1.55, 0, 0.3), 600, cls='mesh thumb', quant=True)
open(OUT+'thumb-flow.svg', 'w').write(svg); print('flow', len(svg))
# ---- apost: porosity on the adapted mesh, around the low-porosity patches ----
svg, _ = svg_from(f['p'], f['t'], f['rho'], 'b', (0.45, 2.25, 0, 0.3), 720, cls='mesh thumb', quant=True)
open(OUT+'thumb-adapt.svg', 'w').write(svg); print('adapt', len(svg))

# ---- poro: Biot + heat, Example 1 (2D) of pp24-02: deformed configuration coloured by p ----
m = MeshTri.init_tensor(np.linspace(0, 1, 17), np.linspace(0, 1, 17))
X, Y = m.p
amp = 1.0
U = 0.1*np.vstack([np.sin(np.pi*X*Y), np.cos(np.pi*X)*np.cos(np.pi*Y)])
pd = m.p + amp*U
c = m.p[:, m.t].mean(1)
pres = np.sin(np.pi*c[0])*np.sin(np.pi*c[1])
outline = '<path class="outline" d="M%.1f %.1fH%.1fV%.1fH%.1fZ"/>'
view = (-0.2, 1.35, -0.35, 1.3)
svg, tr = svg_from(pd, m.t, pres, 'b', view, 480, lo=0, hi=1, cls='mesh thumb')
a = tr(np.array([[0, 1], [0, 1]]))  # corners (0,0) and (1,1)
svg = svg.replace('</svg>', outline % (a[0,0], a[0,1], a[1,0], a[1,1], a[0,0]) + '</svg>')
open(OUT+'thumb-poro.svg', 'w').write(svg); print('poro', len(svg))

# ---- transport: chemotaxis-Navier-Stokes Example 3 of pp22-16: initial cell density ----
m = MeshTri.init_symmetric().scaled([2, 1]).refined(4)
c = m.p[:, m.t].mean(1)
eta0 = sum(70*np.exp(-8*(c[0]-s)**2 - 10*(c[1]-1)**2) for s in (0.2, 0.5, 1.2))
svg, _ = svg_from(m.p, m.t, np.sqrt(eta0), 's', (0, 2, 0, 1), 600, lo=0, hi=np.sqrt(eta0.max()), cls='mesh thumb')
# chemical signal: contour circles of phi0 = 30 exp(-5(x-1)^2 - 5(y-0.5)^2), levels 10, 20
r = lambda lev: np.sqrt(np.log(30/lev)/5)
cx, cy = 300, 150; sc = 300
rings = "".join('<circle class="contour" cx="%d" cy="%d" r="%.1f"/>' % (cx, cy, r(l)*sc) for l in (5, 12, 22))
svg = svg.replace('</svg>', rings + '</svg>')
open(OUT+'thumb-transport.svg', 'w').write(svg); print('transport', len(svg), m.t.shape[1])
