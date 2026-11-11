import numpy as np
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties
T, O, P, I, L = '#33607A', '#E2BE72', '#131C23', '#E1E8EA', '#92B7B4'   # teal, ochre, dark paper, ink, light teal
def svg(body, bg=None, r=12):
    b = f'<rect width="64" height="64" rx="{r}" fill="{bg}"/>' if bg else ''
    return f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 64 64">{b}{body}</svg>'
def text_path(s, size, x, y, fam='STIXGeneral', weight='normal', style='normal'):
    tp = TextPath((0, 0), s, size=size, prop=FontProperties(family=fam, weight=weight, style=style))
    v = tp.vertices; xmin, ymin = v.min(0); xmax, ymax = v.max(0)
    w, h = xmax - xmin, ymax - ymin
    d = []; codes = tp.codes; i = 0
    for (px, py), c in zip(v, codes):
        X = x - w/2 + (px - xmin); Y = y + h/2 - (py - ymin)
        d.append({1: 'M', 2: 'L', 3: 'Q', 4: 'C', 79: 'Z'}[c] + ('' if c == 79 else f'{X:.2f} {Y:.2f}'))
    # collapse repeated Q/C letters handled by SVG implicitly; keep explicit for safety
    out = []; k = 0
    while k < len(d):
        cmd = d[k][0]
        if cmd == 'Q': out.append('Q' + d[k][1:] + ' ' + d[k+1][1:]); k += 2
        elif cmd == 'C': out.append('C' + d[k][1:] + ' ' + d[k+1][1:] + ' ' + d[k+2][1:]); k += 3
        else: out.append(d[k]); k += 1
    return ''.join(out)
opts = {}
# A: red-refined triangle (one level), one child highlighted
A, B, C = (6, 56), (58, 56), (32, 8)
mid = lambda p, q: ((p[0]+q[0])/2, (p[1]+q[1])/2)
ab, bc, ca = mid(A, B), mid(B, C), mid(C, A)
tri = lambda pts, f: '<path d="M%g %gL%g %gL%g %gZ" fill="%s" stroke="%s" stroke-width="3" stroke-linejoin="round"/>' % (*sum(pts, ()), f, P)
opts['A'] = svg(tri((A, ab, ca), T) + tri((ab, B, bc), T) + tri((ca, bc, C), L) + tri((ab, bc, ca), O), P)
# B: Raviart-Thomas element: triangle with outward normal arrows at edge midpoints
def arrow(m, n, L=9):
    x, y = m; ex, ey = x + n[0]*L, y + n[1]*L
    px, py = -n[1], n[0]
    return (f'<path d="M{x:.1f} {y:.1f}L{ex:.1f} {ey:.1f}" stroke="{O}" stroke-width="4" stroke-linecap="round"/>'
            f'<path d="M{ex+n[0]*4:.1f} {ey+n[1]*4:.1f}L{ex-px*4:.1f} {ey-py*4:.1f}L{ex+px*4:.1f} {ey+py*4:.1f}Z" fill="{O}"/>')
A2, B2, C2 = (17, 44), (47, 44), (32, 18)
body = '<path d="M17 44L47 44L32 18Z" fill="%s" stroke="%s" stroke-width="3.5" stroke-linejoin="round"/>' % (T, L)
for p, q in ((A2, B2), (B2, C2), (C2, A2)):
    m = mid(p, q); e = np.array(q) - np.array(p); n = np.array([e[1], -e[0]]); n = n/np.linalg.norm(n)
    cen = (np.array(A2) + np.array(B2) + np.array(C2))/3
    if np.dot(n, np.array(m) - cen) < 0: n = -n
    body += arrow(m, n)
opts['B'] = svg(body, P)
# C: serif monogram CI on teal
opts['C'] = svg(f'<path d="{text_path("CI", 40, 32, 33)}" fill="{I}"/>', T)
# D: monogram over a mesh corner
mesh = ''.join(tri(t, f).replace('stroke-width="3"', 'stroke-width="1.6"') for t, f in [
    (((0, 64), (24, 64), (0, 40)), O), (((24, 64), (0, 40), (22, 42)), '#C79A45'),
    (((24, 64), (22, 42), (44, 64)), '#3D7384'), (((0, 40), (22, 42), (0, 22)), '#3D7384')])
opts['D'] = svg(f'<clipPath id="k"><rect width="64" height="64" rx="12"/></clipPath><g clip-path="url(#k)"><rect width="64" height="64" fill="{T}"/>{mesh}</g>'
                f'<path d="{text_path("CI", 36, 37, 26)}" fill="{I}"/>')
# E: nabla (divergence) in serif
opts['E'] = svg(f'<path d="{text_path("∇", 50, 32, 33, weight="bold")}" fill="{O}"/>', P)
for k, s in opts.items(): pass

