# Camano, Caucao & Inzunza, pseudostress-based fully mixed unsteady Biot-Brinkman (GIANuC2 2026), Example 2:
# layer 70 of SPE10 (Christie & Blunt 2001). Panel (a): log10 of the permeability (exact data).
# Panel (b): steady Darcy velocity magnitude for a unit pressure drop left to right, computed here (TPFA).
import numpy as np, sys, scipy.sparse as sp, scipy.sparse.linalg as spl
sys.path.insert(0, 'scripts'); from sitefig import *
nx, ny, nz = 60, 220, 85
f = open('scripts/SPE10MODEL2_PERM.INC')
vals = []; started = False
for line in f:
    if line.startswith('PERMX'): started = True; continue
    if started:
        if line.startswith('PERMY'): break
        vals.extend(float(v) for v in line.split() if v[0].isdigit() or v[0] == '.')
K = np.array(vals[:nx*ny*nz]).reshape(nz, ny, nx)[69]     # layer 70
print('K range', K.min(), K.max(), 'geo mean', np.exp(np.log(K).mean()))
Lx, Ly = 0.5455, 1.0
hx, hy = Lx/nx, Ly/ny
# transmissibilities (harmonic means), Dirichlet p=1 at x=0, p=0 at x=Lx, no flow top/bottom
idx = lambda i, j: j*nx + i
rows, cols, data = [], [], []; b = np.zeros(nx*ny)
for j in range(ny):
    for i in range(nx):
        k = idx(i, j); diag = 0.0
        for di, dj, h, hh in ((1, 0, hx, hy), (-1, 0, hx, hy), (0, 1, hy, hx), (0, -1, hy, hx)):
            ii, jj = i + di, j + dj
            if 0 <= ii < nx and 0 <= jj < ny:
                Tf = 2*K[j, i]*K[jj, ii]/(K[j, i] + K[jj, ii])*hh/h
                rows.append(k); cols.append(idx(ii, jj)); data.append(-Tf); diag += Tf
            elif dj == 0:
                Tf = 2*K[j, i]*hh/h; diag += Tf
                if ii < 0: b[k] += Tf*1.0
        rows.append(k); cols.append(k); data.append(diag)
A = sp.csr_matrix((data, (rows, cols)), shape=(nx*ny, nx*ny))
p = spl.spsolve(A, b).reshape(ny, nx)
# cell velocity from face fluxes
ux = np.zeros((ny, nx + 1)); uy = np.zeros((ny + 1, nx))
for j in range(ny):
    for i in range(1, nx):
        ux[j, i] = -2*K[j, i]*K[j, i-1]/(K[j, i] + K[j, i-1])*(p[j, i] - p[j, i-1])/hx
    ux[j, 0] = -2*K[j, 0]*(p[j, 0] - 1.0)/hx; ux[j, nx] = -2*K[j, -1]*(0.0 - p[j, -1])/hx
for j in range(1, ny):
    uy[j] = -2*K[j]*K[j-1]/(K[j] + K[j-1])*(p[j] - p[j-1])/hy
vel = np.hypot(0.5*(ux[:, 1:] + ux[:, :-1]), 0.5*(uy[1:] + uy[:-1]))
xe = np.linspace(0, Lx, nx + 1); ye = np.linspace(0, Ly, ny + 1)   # j=0 at y=0 (bottom)
W = 190; gap = 34
P1 = Panel((0, Lx, 0, Ly), W, 0, 0); P2 = Panel((0, Lx, 0, Ly), W, W + gap, 0)
body = grid_runs(P1, bands(np.log10(K), quant=True), xe, ye, 'b')
body += grid_runs(P2, bands(np.log10(vel + 1e-30), quant=True), xe, ye, 's')
H = P1.h
body += label(W/2, H + 22, 'log₁₀ K') + label(W + gap + W/2, H + 22, 'log₁₀ |velocity|')
standalone('assets/img/figs/spe10-layer70.svg', 2*W + gap, H + 30, body, grid=True)
print('svg written')
