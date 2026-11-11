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

# ---------- thumbnails ----------
def thumb(name, P, fieldf, lo, hi, inside=None, transform=None):
    w, h = 240, 160
    T = Delaunay(P).simplices
    c = P[T].mean(1)
    Q = P if transform is None else transform(P)
    s = '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" class="mesh thumb" aria-hidden="true" focusable="false">' % (w, h)
    s += paths(Q, T, bands(fieldf(c), lo, hi), inside) + '</svg>'
    open('_includes/thumb-%s.svg' % name, 'w').write(s)
    print(name, len(T), len(s))

w, h = 240, 160
# poroelasticity: block loaded on top, deformed, pressure field
xs, ys = np.meshgrid(np.linspace(14, 226, 15), np.linspace(10, 132, 9))
P = np.c_[xs.ravel(), ys.ravel()] + rng.uniform(-3, 3, (xs.size, 2))*0
def deform(P):
    x, y = (P[:,0]-14)/212, (P[:,1]-10)/122
    Q = P.copy(); Q[:,1] += 18*np.sin(np.pi*x)*(1-y)**1.2; Q[:,0] += 6*(x-0.5)*(1-y)
    return Q
thumb('poro', P, lambda c: np.sin(np.pi*(c[:,0]-14)/212)*(1-(c[:,1]-10)/122), 0, 1, transform=deform)
# flow around inclusion
P = jitter_grid(w, h, 15)
r = lambda p: np.hypot(p[:,0]-110, p[:,1]-80)
P = P[r(P) > 30]
th = np.linspace(0, 2*np.pi, 28, endpoint=False)
P = np.vstack([P, np.c_[110+30*np.cos(th), 80+30*np.sin(th)], np.c_[110+38*np.cos(th+0.1), 80+38*np.sin(th+0.1)]])
def speed(c):
    x, y = c[:,0]-110, c[:,1]-80; rr = np.hypot(x, y)+1e-9; a = 30
    ur = (1-a**2/rr**2)*x/rr; ut = -(1+a**2/rr**2)*y/rr
    return np.hypot(ur, ut)
thumb('flow', P, speed, 0.2, 1.9, inside=lambda c: np.hypot(c[0]-110, c[1]-80) > 30.5)
# transport plume
P = jitter_grid(w, h, 14)
def plume(c):
    x, y = c[:,0]/w, c[:,1]/h
    sig = 0.04+0.25*x
    return np.exp(-(y-0.5-0.12*np.sin(5*x))**2/(2*sig**2)) * np.exp(-1.2*x) + 0.0
thumb('transport', P, plume, 0.02, 0.95)
# a posteriori: L-shaped domain refined toward reentrant corner
pts = [(20,20),(220,20),(220,80),(120,80),(120,140),(20,140)]
corner = np.array([120, 80])
P = jitter_grid(w, h, 26, 0.25)
for hh, rad in [(13, 50), (6, 22), (3, 9)]:
    G = jitter_grid(w, h, hh, 0.25); P = np.vstack([P, G[np.hypot(*(G-corner).T) < rad]])
inL = lambda p: (20 <= p[0] <= 220 and 20 <= p[1] <= 140) and not (p[0] > 120 and p[1] > 80)
bd = []
for (a, b) in zip(pts, pts[1:]+pts[:1]):
    n = int(np.hypot(b[0]-a[0], b[1]-a[1])/8)
    for t in np.linspace(0, 1, n, endpoint=False): bd.append((a[0]+t*(b[0]-a[0]), a[1]+t*(b[1]-a[1])))
bd = np.array(bd)
P = np.array([p for p in P if inL(p) and min(abs(p[0]-20),abs(p[0]-220),abs(p[1]-20),abs(p[1]-140)) > 3 and not (abs(p[0]-120)<3 and p[1]>80) and not (abs(p[1]-80)<3 and p[0]>120)])
P = np.vstack([P, bd])
thumb('adapt', P, lambda c: (np.hypot(*(c-corner).T)+1)**(-1/3), 0.12, 0.5, inside=inL)
