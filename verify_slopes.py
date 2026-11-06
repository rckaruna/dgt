"""verify_slopes.py -- independent numpy port of .dstudy_slopes_draws() (dgt 0.7.0).

Checks (all must PASS):
  1. link-scale curve equals 1 - V_kk(m)/sigma_k^2 with V = (m Ibar + Sigma^-1)^-1
  2. both curves are 0 at m = 0, non-decreasing in m, and -> 1 as m -> inf
  3. response-scale <= link-scale at every m (logit)
  4. with one slope, S = 1, intercept only: Ibar = E[p(1-p)] at the intercept,
     so link curve = s^2 / (s^2 + 1/(m w)) -- the random-intercept closed form
     with resid variance 1/w in place of pi^2/3
Prints reference values that tests/testthat/test-dstudy-slopes.R compares to.
"""
import numpy as np

def gh_nodes(n=40):
    k = np.arange(1, n); J = np.zeros((n, n))
    J[k - 1, k] = np.sqrt(k); J[k, k - 1] = np.sqrt(k)
    ev, V = np.linalg.eigh(J); o = np.argsort(ev)
    return ev[o], V[0, o] ** 2

def plogis(x): return 1 / (1 + np.exp(-x))

def slopes_core(beta, sd, X, S, Xp, Sp, n_grid, nodes=40):
    z, w = gh_nodes(nodes); k = S.shape[1]
    p = plogis(X @ beta); wt = p * (1 - p)
    Ibar = (S * wt[:, None]).T @ S / S.shape[0]
    Sinv = np.diag(1 / sd ** 2)
    a = Xp @ beta; cc = np.diag(Sp)
    var_pi = np.array([(w * plogis(a[j] + cc[j] * sd[j] * z) ** 2).sum() - (w * plogis(a[j] + cc[j] * sd[j] * z)).sum() ** 2 for j in range(k)])
    link = np.zeros((len(n_grid), k)); resp = np.zeros((len(n_grid), k))
    for i, m in enumerate(n_grid):
        V = np.diag(sd ** 2) if m == 0 else np.linalg.inv(m * Ibar + Sinv)
        v = np.clip(np.diag(V), 0, sd ** 2)
        link[i] = 1 - v / sd ** 2
        for j in range(k):
            tau = np.sqrt(sd[j] ** 2 - v[j])
            if tau <= 0: resp[i, j] = 0; continue
            th = tau * z
            Epi = np.array([(w * plogis(a[j] + cc[j] * (t + np.sqrt(v[j]) * z))).sum() for t in th])
            resp[i, j] = ((w * Epi ** 2).sum() - (w * Epi).sum() ** 2) / var_pi[j]
    return link, resp, Ibar

ok = True
def check(name, cond):
    global ok; ok &= bool(cond); print(("PASS " if cond else "FAIL ") + name)

# ---- reference case: 2 slopes, 6 tasks, intercept + 2 fixed covariates ----
X  = np.array([[1, 1, 0], [1, 0, 1], [1, -1, 0], [1, 0, -1], [1, 1, 1], [1, -1, 1]], float)
S  = X[:, 1:].copy()                      # slopes on the same covariates
beta = np.array([0.3, -0.8, -1.2]); sd = np.array([0.9, 1.4])
Xp = np.array([[1, 1, 0], [1, 0, 1]], float); Sp = np.array([[1, 0], [0, 1]], float)
n_grid = np.array([0, 1, 2, 5, 10, 20, 50, 200, 5000])
link, resp, Ibar = slopes_core(beta, sd, X, S, Xp, Sp, n_grid)

# 1. closed form
p = plogis(X @ beta); wt = p * (1 - p)
I2 = (S * wt[:, None]).T @ S / 6
for i, m in enumerate(n_grid):
    V = np.diag(sd**2) if m == 0 else np.linalg.inv(m * I2 + np.diag(1 / sd**2))
    assert np.allclose(link[i], 1 - np.diag(V) / sd**2, atol=1e-12)
check("1 link curve = 1 - V_kk/sigma^2", True)
# 2. limits and monotonicity
check("2a both curves 0 at m=0", np.allclose(link[0], 0) and np.allclose(resp[0], 0))
check("2b non-decreasing in m", (np.diff(link, axis=0) >= -1e-12).all() and (np.diff(resp, axis=0) >= -1e-9).all())
check("2c -> 1 at m=5000", (link[-1] > 0.995).all() and (resp[-1] > 0.99).all())
# 3. ordering
check("3 response <= link at every m", (resp <= link + 1e-9).all())
# 4. one slope, S=1, reduces to the intercept closed form
X1 = np.array([[1.0]]); S1 = np.array([[1.0]]); b1 = np.array([0.4]); s1 = np.array([1.1])
l1, r1, _ = slopes_core(b1, s1, X1, S1, np.array([[1.0]]), np.array([[1.0]]), n_grid)
w0 = plogis(0.4) * (1 - plogis(0.4))
exp1 = np.array([0 if m == 0 else s1[0]**2 / (s1[0]**2 + 1 / (m * w0)) for m in n_grid])
check("4 single slope matches intercept closed form", np.allclose(l1[:, 0], exp1, atol=1e-12))

print("\nReference values for test-dstudy-slopes.R (reference case, m = 20):")
print("link     =", np.round(link[n_grid == 20][0], 10).tolist())
print("response =", np.round(resp[n_grid == 20][0], 10).tolist())
print("Ibar     =", np.round(Ibar, 10).tolist())
print("\nALL PASS" if ok else "\nSOME CHECKS FAILED")
