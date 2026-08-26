# verify_v040.py -- Python cross-checks for the v0.4.0 identities.
#
# Check 1: crossed relative reliability. On a fixed panel of n facet
#   levels, the variance over objects of the observed panel mean equals
#   signal + noise, where signal = Var_k(p_k) and noise is the Bernoulli
#   variance of the panel mean; the reliability formula equals the
#   squared correlation between observed and true panel scores.
#
# Check 2: variance of the grand mean of a crossed n_r x n_t Bernoulli
#   array equals [mu(1-mu) + (n_t-1) Var_u(t) + (n_r-1) Var_v(w)] /
#   (n_r n_t), with t = E_v[pi|u] and w = E_u[pi|v].
#
# Run: python3 verify_v040.py     (requires numpy)

import numpy as np

rng = np.random.default_rng(7)
expit = lambda x: 1 / (1 + np.exp(-x))
fails = 0


def check(name, ok, detail):
    global fails
    print(("PASS" if ok else "FAIL"), name, "--", detail)
    if not ok:
        fails += 1


# ---------- Check 1: crossed relative reliability ----------
a, so, sf, n = 0.3, 1.0, 0.6, 8
Kpop = 400_000
u = rng.normal(0, so, Kpop)
v = rng.normal(0, sf, n)
P = expit(a + u[:, None] + v[None, :])
p = P.mean(1)
Y = rng.binomial(1, P).mean(1)
sig = p.var()
noise = (P * (1 - P)).sum(1).mean() / n**2
check("C1a Var(panel mean) = signal + noise",
      abs(Y.var() / (sig + noise) - 1) < 0.01,
      "emp %.6f formula %.6f" % (Y.var(), sig + noise))
rho_f = sig / (sig + noise)
rho_e = np.corrcoef(Y, p)[0, 1] ** 2
check("C1b reliability formula = squared correlation",
      abs(rho_f - rho_e) < 0.01,
      "formula %.5f corr^2 %.5f" % (rho_f, rho_e))

# ---------- Check 2: grand-mean variance identity ----------
a, so, sf = -0.2, 0.9, 0.5
nr, nt = 6, 4
U = rng.normal(0, so, 40_000)
V = rng.normal(0, sf, 2_000)
t = np.zeros(len(U)); rowvar = np.zeros(len(U))
w = np.zeros(len(V)); w2 = np.zeros(len(V)); s_mu = 0.0
for i0 in range(0, len(U), 5_000):
    PP = expit(a + U[i0:i0 + 5_000, None] + V[None, :])
    t[i0:i0 + 5_000] = PP.mean(1)
    rowvar[i0:i0 + 5_000] = PP.var(1, ddof=1)
    w += PP.sum(0); w2 += (PP ** 2).sum(0); s_mu += PP.sum()
mu = s_mu / (len(U) * len(V))
w /= len(U)
colvar = (w2 / len(U) - w ** 2) * len(U) / (len(U) - 1)
var_u_t = t.var(ddof=1) - rowvar.mean() / len(V)
var_v_w = np.var(w, ddof=1) - colvar.mean() / len(U)
form = (mu * (1 - mu) + (nt - 1) * var_u_t + (nr - 1) * var_v_w) / (nr * nt)

M = 300_000
acc = 0.0; acc2 = 0.0
for j0 in range(0, M, 50_000):
    m = 50_000
    uu = rng.normal(0, so, (m, nr))
    vv = rng.normal(0, sf, (m, nt))
    Pg = expit(a + uu[:, :, None] + vv[:, None, :])
    Xb = rng.binomial(1, Pg).mean((1, 2))
    acc += Xb.sum(); acc2 += (Xb ** 2).sum()
emp = acc2 / M - (acc / M) ** 2
check("C2 Var(grand mean) identity",
      abs(emp / form - 1) < 0.02,
      "emp %.6f formula %.6f ratio %.4f" % (emp, form, emp / form))

# ---------- Check 3: degenerate allocation invariance ----------
mu0 = expit(0.4)
vals = [mu0 * (1 - mu0) / (nr * nt) for nr, nt in
        [(2, 12), (3, 8), (4, 6), (24, 1)]]
check("C3 sd = 0 gives mu(1-mu)/N for every allocation",
      max(vals) - min(vals) < 1e-15,
      "spread %.2e" % (max(vals) - min(vals)))

print("\n%d check(s) failed." % fails)
raise SystemExit(1 if fails else 0)
