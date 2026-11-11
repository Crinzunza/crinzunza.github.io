# Camano, Caucao & Inzunza, five-field fully dynamic Biot-Brinkman (M2AN, accepted), Example 1:
# exact solution at t = 0 on the unit square: sigma_p(1,2), |u|, |eta| and p.
import numpy as np, sys
sys.path.insert(0, 'scripts'); from sitefig import *
from skfem import MeshTri
m = MeshTri.init_tensor(np.linspace(0, 1, 25), np.linspace(0, 1, 25))
x, y = m.p[:, m.t].mean(1); pi = np.pi
s12 = -2*pi*np.sin(pi*x)*np.sin(pi*y)                          # 2 mu e_12(eta), mu = 1 scaling
umag = np.hypot(np.sin(pi*x)*np.sin(2*pi*y), np.sin(2*pi*x)*np.sin(pi*y))
emag = np.hypot(np.sin(pi*x)*np.cos(pi*y), np.cos(pi*x)*np.sin(pi*y))
p = np.sin(pi*x)*np.cos(pi*y/2) - 2/pi
W, gap = 150, 22; body = ''
for k, (v, pal, sym, lab) in enumerate([(s12, 'b', True, 'σ₁₂'), (umag, 's', False, '|u|'), (emag, 's', False, '|η|'), (p, 'b', True, 'p')]):
    P = Panel((0, 1, 0, 1), W, k*(W + gap), 0)
    body += tri_paths(P, m.p, m.t, bands(v, sym=sym) if sym else bands(v, lo=0), pal)
    body += label(k*(W + gap) + W/2, W + 22, lab)
standalone('assets/img/figs/fivefield-exact.svg', 4*W + 3*gap, W + 30, body)
