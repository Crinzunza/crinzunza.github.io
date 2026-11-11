import numpy as np, sys
sys.argv = ['x']; exec(open('scripts/favicon_base.py').read().split("opts = {}")[0])
def arr(m, n, L, w=4, head=4.5, col=O):
    x, y = m; ex, ey = x + n[0]*L, y + n[1]*L; px, py = -n[1], n[0]
    return (f'<path d="M{x:.1f} {y:.1f}L{ex:.1f} {ey:.1f}" stroke="{col}" stroke-width="{w}" stroke-linecap="round"/>'
            f'<path d="M{ex+n[0]*head:.1f} {ey+n[1]*head:.1f}L{ex-px*head:.1f} {ey-py*head:.1f}L{ex+px*head:.1f} {ey+py*head:.1f}Z" fill="{col}"/>')
def tri_rt(A, B, C, L):
    out = ''; cen = (np.array(A) + np.array(B) + np.array(C))/3
    for p, q in ((A, B), (B, C), (C, A)):
        m = (np.array(p) + np.array(q))/2; e = np.array(q) - np.array(p); n = np.array([e[1], -e[0]], float); n /= np.linalg.norm(n)
        if np.dot(n, m - cen) < 0: n = -n
        out += arr(m, n, L)
    return out
v = {}
# F1: dark background, teal triangle with CI inside, ochre outward fluxes
A, B, C = (10, 50), (54, 50), (32, 11)
v['F1'] = svg(f'<path d="M10 50L54 50L32 11Z" fill="{T}" stroke="{L}" stroke-width="2.5" stroke-linejoin="round"/>'
              + tri_rt(A, B, C, 6) + f'<path d="{text_path("CI", 19, 32, 39)}" fill="{I}"/>', P)
# F2: the square tile is the element: CI with four outward normal fluxes (RT on a quadrilateral)
body = f'<rect x="10" y="10" width="44" height="44" rx="3" fill="{T}" stroke="{L}" stroke-width="2.5"/>'
for m, n in (((32, 10), (0, -1)), ((54, 32), (1, 0)), ((32, 54), (0, 1)), ((10, 32), (-1, 0))):
    body += arr(m, n, 3.5, w=3.5, head=4)
body += f'<path d="{text_path("CI", 27, 32, 32.5)}" fill="{I}"/>'
v['F2'] = svg(body, P)
# F3: monogram tile (C) with a small RT element as the dot over the I
body = f'<path d="{text_path("CI", 36, 30, 37)}" fill="{I}"/>'
A3, B3, C3 = (43, 22), (57, 22), (50, 10)
body += f'<path d="M43 22L57 22L50 10Z" fill="{O}"/>'
v['F3'] = svg(body, T)
for k, s in v.items(): open('assets/img/favicon.svg' if k == 'F1' else '/dev/null', 'w').write(s)
