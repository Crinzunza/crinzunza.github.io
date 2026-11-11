"""Helpers that turn meshes and fields into site-styled SVG (bands coloured by CSS variables)."""
import numpy as np
NB = 7

def bands(val, lo=None, hi=None, quant=False, sym=False):
    val = np.asarray(val, float)
    if quant:
        edges = np.quantile(val, np.linspace(0, 1, NB + 1)[1:-1]); return np.searchsorted(edges, val)
    lo = np.nanmin(val) if lo is None else lo; hi = np.nanmax(val) if hi is None else hi
    if sym: m = max(abs(lo), abs(hi)); lo, hi = -m, m
    return np.clip(((val - lo) / (hi - lo) * NB).astype(int), 0, NB - 1)

class Panel:
    """Maps model coordinates in `view`=(x0,x1,y0,y1) to an SVG box at offset (ox, oy) of width w."""
    def __init__(self, view, w, ox=0, oy=0):
        self.x0, self.x1, self.y0, self.y1 = view; self.s = w / (self.x1 - self.x0)
        self.w, self.h, self.ox, self.oy = w, (self.y1 - self.y0) * self.s, ox, oy
    def tr(self, P):
        P = np.asarray(P, float)
        return np.c_[(P[0] - self.x0) * self.s + self.ox, (self.y1 - P[1]) * self.s + self.oy]

def fmt(v):
    s = "%.1f" % v
    return s[:-2] if s.endswith('.0') else s

def tri_paths(panel, p, t, q, pal, skip=None):
    P = panel.tr(p); out = {}
    for k in range(t.shape[1]):
        if skip is not None and skip[k]: continue
        a = P[t[:, k]]
        out.setdefault(int(q[k]), []).append("M%sL%sL%sZ" % tuple(fmt(x) + " " + fmt(y) for x, y in a))
    return "".join('<path class="%s%d" d="%s"/>' % (pal, b, "".join(v)) for b, v in sorted(out.items()))

def grid_runs(panel, q, x_edges, y_edges, pal, outline=False):
    """q: (ny, nx) band indices on a tensor grid; merges runs horizontally, then identical runs vertically."""
    open_rects = {}; done = []
    for j in range(q.shape[0]):
        runs = []; i = 0
        while i < q.shape[1]:
            k = i
            while k + 1 < q.shape[1] and q[j, k + 1] == q[j, i]: k += 1
            runs.append((i, k + 1, int(q[j, i]))); i = k + 1
        new_open = {}
        for r in runs:
            if r in open_rects: new_open[r] = (open_rects.pop(r)[0], j + 1)
            else: new_open[r] = (j, j + 1)
        done += [(r, span) for r, span in open_rects.items()]; open_rects = new_open
    done += list(open_rects.items())
    out = {}
    for (i0, i1, b), (j0, j1) in done:
        (xa, ya), (xb, yb) = panel.tr([[x_edges[i0], x_edges[i1]], [y_edges[j1], y_edges[j0]]])
        out.setdefault(b, []).append("M%s %sH%sV%sH%sZ" % (fmt(xa), fmt(ya), fmt(xb), fmt(yb), fmt(xa)))
    return "".join('<path class="%s%d" d="%s"/>' % (pal, b, "".join(v)) for b, v in sorted(out.items()))

def svg(width, height, body, cls="mesh thumb"):
    return ('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" class="%s" aria-hidden="true" focusable="false">'
            % (round(width), round(height), cls)) + body + '</svg>'

def label(x, y, text, anchor="middle"):
    return '<text class="panel-label" x="%s" y="%s" text-anchor="%s">%s</text>' % (fmt(x), fmt(y), anchor, text)

STYLE = """<style>
svg{--b0:#E2BE72;--b1:#EAD7A6;--b2:#E6E5D7;--b3:#C5D9D3;--b4:#92B7B4;--b5:#5E8C99;--b6:#33607A;
--s0:#E9EEEA;--s1:#D2E0DB;--s2:#B2CDC8;--s3:#8DB4B4;--s4:#6A9BA6;--s5:#44788F;--s6:#24506A;--line:#F6F8F7;--ink:#56666E}
@media (prefers-color-scheme: dark){svg{--b0:#8C6A2A;--b1:#6B5C3C;--b2:#384447;--b3:#2D474C;--b4:#2E5862;--b5:#3D7384;--b6:#5C95A8;
--s0:#1E2B33;--s1:#23363F;--s2:#2A4550;--s3:#315A66;--s4:#3B7182;--s5:#5290A3;--s6:#86BED0;--line:#131C23;--ink:#9BABB1}}
path{stroke:var(--line);stroke-width:.3;stroke-linejoin:round}
.grid path{stroke:none;shape-rendering:crispEdges}
.b0{fill:var(--b0)}.b1{fill:var(--b1)}.b2{fill:var(--b2)}.b3{fill:var(--b3)}.b4{fill:var(--b4)}.b5{fill:var(--b5)}.b6{fill:var(--b6)}
.s0{fill:var(--s0)}.s1{fill:var(--s1)}.s2{fill:var(--s2)}.s3{fill:var(--s3)}.s4{fill:var(--s4)}.s5{fill:var(--s5)}.s6{fill:var(--s6)}
text{font-family:'STIX Two Text','Times New Roman',serif;font-size:15px;fill:var(--ink)}
.contour,.outline{fill:none;stroke:var(--ink);stroke-width:1;opacity:.5}.outline{stroke-dasharray:4 3}
</style>"""

def standalone(path, width, height, body, grid=False):
    open(path, 'w').write('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d"%s>' % (round(width), round(height), ' class="grid"' if grid else '') + STYLE + body + '</svg>')
