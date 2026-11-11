# Regenerates the mesh figures. Run from the repository root: python3 scripts/generate_meshes.py
import numpy as np
from scipy.spatial import Delaunay
rng = np.random.default_rng(7)
NB = 7

def bands(vals, lo, hi):
    q = np.clip(((vals-lo)/(hi-lo)*NB).astype(int), 0, NB-1)
    return q

def paths(pts, tri, q, inside=None):
    out = {}
    for t, b in zip(tri, q):
        if inside is not None and not inside(pts[t].mean(0)): continue
        p = pts[t]
        d = "M%d %dL%d %dL%d %dZ" % tuple(np.round(p.ravel()).astype(int))
        out.setdefault(int(b), []).append(d)
    return "".join('<path class="b%d" d="%s"/>' % (b, "".join(v)) for b, v in sorted(out.items()))

def jitter_grid(W, H, h, j=0.35):
    xs = np.arange(-h, W+2*h, h); ys = np.arange(-h, H+2*h, h*0.87)
    P = []
    for i, y in enumerate(ys):
        for x in xs:
            P.append((x + (h/2 if i % 2 else 0) + rng.uniform(-j, j)*h, y + rng.uniform(-j, j)*h))
    return np.array(P)

def near(P, f, d):
    return P[np.abs(P[:,1]-f(P[:,0])) < d]

# ---------- hero ----------
W, H = 1440, 440
curve = lambda x: 0.52*H + 55*np.sin(2*np.pi*x/720 + 0.6) + 18*np.sin(2*np.pi*x/260)
def field(c):
    x, y = c[:,0]/W, c[:,1]/H
    s = curve(c[:,0])
    above = (c[:,1] < s)
    v = 0.9*np.exp(-((x-0.78)**2/0.02 + (y-0.3)**2/0.08)) - 0.8*np.exp(-((x-0.22)**2/0.03 + (y-0.75)**2/0.06))
    v += 0.25*np.sin(3*np.pi*x)*np.cos(2*np.pi*y) + np.where(above, 0.12, -0.12)
    return v
levels = []
P0 = jitter_grid(W, H, 92)
P1 = np.vstack([P0, near(jitter_grid(W, H, 40), curve, 95)])
P2 = np.vstack([P1, near(jitter_grid(W, H, 17), curve, 34)])
svg = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" preserveAspectRatio="xMidYMid slice" class="mesh" aria-hidden="true" focusable="false">' % (W, H)]
for k, P in enumerate([P0, P1, P2]):
    T = Delaunay(P).simplices
    c = P[T].mean(1)
    print("hero level", k, len(T))
    svg.append('<g class="lvl lvl%d">%s</g>' % (k, paths(P, T, bands(field(c), -0.85, 0.95))))
svg.append('</svg>')
open('_includes/hero-mesh.svg', 'w').write("".join(svg))

