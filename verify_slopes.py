"""verify_slopes.py -- independent numpy port of .dstudy_slopes_draws() (dgt 0.7.1).

Checks (all must PASS):
  1. link-scale curve equals 1 - V_kk(m)/sigma_k^2 with V = (m Ibar + Sigma^-1)^-1
  2. both curves are 0 at m = 0, non-decreasing in m, and -> 1 as m -> inf
  3. response-scale <= link-scale at every m (logit)
  4. with one slope, S = 1, intercept only: Ibar = E[p(1-p)] at the intercept,
     so link curve = s^2 / (s^2 + 1/(m w)) -- the random-intercept closed form
     with resid variance 1/w in place of pi^2/3
  5. method = "simulate": 0 at m = 0, link monotone, below the plug-in curve
  6. against respondents simulated from the model and estimated by posterior
     mean, "simulate" matches the squared true/estimated correlation and
     "plugin" overstates it
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

# ---- method = "simulate" (v0.7.1): average over respondents and task draws ----
def slopes_core_sim(beta, sd, X, S, Xp, Sp, n_grid, K=2000, seed=1, nodes=40):
    rng = np.random.default_rng(seed); z, w = gh_nodes(nodes); k = S.shape[1]; nR = len(S)
    fx = X @ beta; Sinv = np.diag(1 / sd ** 2); a = Xp @ beta; cc = np.diag(Sp); M = int(max(n_grid))
    theta = rng.normal(size=(K, k)) * sd; v = np.zeros((K, len(n_grid), k))
    pos = {int(m): i for i, m in enumerate(n_grid)}
    for r in range(K):
        if 0 in pos: v[r, pos[0]] = sd ** 2
        order = np.concatenate([rng.permutation(nR) for _ in range(-(-M // nR))])[:M]
        p = plogis(fx[order] + S[order] @ theta[r]); wr = p * (1 - p); info = Sinv.copy()
        for t in range(M):
            s_ = S[order[t]]; info += wr[t] * np.outer(s_, s_)
            if t + 1 in pos: v[r, pos[t + 1]] = np.minimum(np.diag(np.linalg.inv(info)), sd ** 2)
    var_pi = np.array([(w * plogis(a[j] + cc[j] * sd[j] * z) ** 2).sum() - (w * plogis(a[j] + cc[j] * sd[j] * z)).sum() ** 2 for j in range(k)])
    link = np.clip(1 - v.mean(0) / sd ** 2, 0, 1); resp = np.zeros_like(link)
    for i in range(len(n_grid)):
        for j in range(k):
            g = plogis(a[j] + cc[j] * (np.sqrt(link[i, j]) * theta[:, j, None] + np.sqrt(v[:, i, j, None]) * z[None, :]))
            resp[i, j] = np.clip(1 - np.mean((g ** 2) @ w - (g @ w) ** 2) / var_pi[j], 0, 1)
    return link, resp

def simulate_truth(beta, sd, X, S, Xp, Sp, m, R=1000, nis=4000, seed=3):
    """Respondents simulated from the model; each one's coefficients and probe
    probability estimated by posterior mean, with the population parameters
    known. Importance sampling with the prior as proposal: slower than a
    Laplace proposal but it cannot miss a skewed or heavy-tailed posterior,
    which a mode-centred proposal does when choices are near-deterministic.
    Returns squared correlations between true and estimated theta_k (link) and
    between true and estimated probe probability (response), and the median
    effective sample size."""
    rng = np.random.default_rng(seed); k = S.shape[1]; fx = X @ beta
    a = Xp @ beta; cc = np.diag(Sp); th = rng.normal(size=(R, k)) * sd
    eh = np.zeros((R, k)); pt = plogis(a + cc * th); ph = np.zeros((R, k)); ess = np.zeros(R)
    for i in range(R):
        idx = rng.choice(len(X), m, replace=(m > len(X))); Si = S[idx]; f = fx[idx]
        y = (rng.random(m) < plogis(f + Si @ th[i])).astype(float)
        d = rng.normal(size=(nis, k)) * sd; lin = f[None, :] + d @ Si.T
        ll = (y * lin - np.logaddexp(0, lin)).sum(1); wt = np.exp(ll - ll.max()); wt /= wt.sum()
        eh[i] = wt @ d; ph[i] = wt @ plogis(a + cc * d); ess[i] = 1 / (wt ** 2).sum()
    r2 = lambda u, v_: np.array([np.corrcoef(u[:, j], v_[:, j])[0, 1] ** 2 for j in range(k)])
    return r2(th, eh), r2(pt, ph), np.median(ess)

ls, rs = slopes_core_sim(beta, sd, X, S, Xp, Sp, n_grid, K=3000, seed=1)
check("5a simulate: both curves 0 at m=0", np.allclose(ls[0], 0) and np.allclose(rs[0], 0, atol=1e-10))
check("5b simulate: link non-decreasing in m, -> 1", (np.diff(ls, axis=0) >= -1e-12).all() and (ls[-1] > 0.99).all() and (rs[-1] > 0.98).all())
check("5c simulate link < plug-in link for m >= 1 (Jensen)", (ls[1:] < link[1:]).all())
big = np.array([2.5, 3.0]); lb, _, _ = slopes_core(beta, big, X, S, Xp, Sp, np.array([20])); lsb, _ = slopes_core_sim(beta, big, X, S, Xp, Sp, np.array([20]), K=3000, seed=2)
gap_ref = (link[n_grid == 20][0] - ls[n_grid == 20][0]).mean(); gap_big = (lb[0] - lsb[0]).mean()
check("5d the plug-in overstatement grows with heterogeneity (%.3f -> %.3f)" % (gap_ref, gap_big), gap_big > gap_ref)

# 6. against respondents simulated from the model: a paired design with 5 slopes,
#    in a harsh regime (strong effects, many near-deterministic choices) and a mild one
rng = np.random.default_rng(10); nset = 40; kk = 5
Sd = rng.choice([-1.0, 0.0, 1.0], size=(nset, kk), p=[0.4, 0.2, 0.4]) * rng.uniform(0.6, 1.6, size=(nset, kk))
for scale, label in [(1.0, "harsh"), (0.4, "mild")]:
    bd = -np.ones(kk) * scale; sdd = np.array([1.0, 1.2, 0.8, 1.5, 1.3]) * scale
    Xpd = np.eye(kk) * 1.4; Spd = np.eye(kk) * 1.4
    lp, rp, _ = slopes_core(bd, sdd, Sd, Sd, Xpd, Spd, np.array([20]))
    lsim, rsim = slopes_core_sim(bd, sdd, Sd, Sd, Xpd, Spd, np.array([20]), K=3000, seed=4)
    lt, rt, ess = simulate_truth(bd, sdd, Sd, Sd, Xpd, Spd, 20)
    print("   %s, m=20  link: plug-in %.3f  simulate %.3f  simulated respondents %.3f | response: plug-in %.3f  simulate %.3f  simulated respondents %.3f | median ESS %d"
          % (label, lp.mean(), lsim.mean(), lt.mean(), rp.mean(), rsim.mean(), rt.mean(), ess))
    check("6 %s: simulate link within 0.06 per slope and 0.025 on average" % label, np.abs(lsim[0] - lt).max() < 0.06 and abs(lsim.mean() - lt.mean()) < 0.025)
    check("6 %s: simulate response within 0.04 on average" % label, abs(rsim.mean() - rt.mean()) < 0.04)
    check("6 %s: plug-in link is further from the simulated respondents than simulate is" % label, abs(lp.mean() - lt.mean()) > abs(lsim.mean() - lt.mean()))

lref, rref = slopes_core_sim(beta, sd, X, S, Xp, Sp, np.array([20]), K=20000, seed=9)
print("\nReference values for test-dstudy-slopes.R, method = simulate (reference case, m = 20, K = 20000):")
print("SIM link     =", np.round(lref[0], 4).tolist())
print("SIM response =", np.round(rref[0], 4).tolist())

print("\nReference values for test-dstudy-slopes.R, method = plugin (reference case, m = 20):")
print("link     =", np.round(link[n_grid == 20][0], 10).tolist())
print("response =", np.round(resp[n_grid == 20][0], 10).tolist())
print("Ibar     =", np.round(Ibar, 10).tolist())
print("\nALL PASS" if ok else "\nSOME CHECKS FAILED")
