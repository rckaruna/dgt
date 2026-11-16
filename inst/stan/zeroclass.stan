// Zero-class binomial GLMM (dgt). Written for flaky tests: objects are tests,
// groups are projects, conditions are machine configurations.
//
// A test either never fails (class 0) or is flaky-capable (class 1).
// P(class 1) = inv_logit(gamma + a_p) for a test of project p.
// In class 1, failures in cell (test t, configuration c) are
//   y_tc ~ Binomial(n_tc, inv_logit(alpha + b_p + u_t + v_c + w_pc + e_tc)),
// with u_t ~ N(0, sd_u^2), v_c ~ N(0, sd_v^2), w_pc ~ N(0, sd_w^2),
// e_tc ~ N(0, sd_e^2), and (a_p, b_p) bivariate normal with correlation rho.
//
// Tests with at least one failure are in class 1; their linear predictor
// theta_t = alpha + b_p + u_t is sampled (hierarchically centred) and
// e_tc is integrated out: by adaptive Gauss-Hermite in cells with failures,
// through the log J table in cells without.
// Tests with no failure are exchangeable within a project, so they enter as a
// count K_p times the log probability that one such test shows no failure:
//   log( 1 - pi_p + pi_p * I_p ),
//   I_p = int prod_c J(alpha + b_p + v_c + w_pc + x) N(x | 0, sd_u^2) dx,
//   J(eta) = E_e[(1 - inv_logit(eta + e))^n].
// log J is tabulated on a grid at every evaluation (windowed Gauss-Legendre
// with a closed-form left tail) and read off by cubic Hermite interpolation;
// the outer integral uses windowed Gauss-Legendre with a closed-form left tail.
functions {
  // log J(t) = log E_z[(1 - inv_logit(t + se z))^n], z ~ N(0, 1), and its
  // derivative in t. Gauss-Legendre on the window where the expected failure
  // count n p runs from exp(-12) to exp(6), clipped to |z| <= 8 unless the
  // whole window lies below -8; first-order closed form for the left tail.
  vector logJ2(real t, real n, real se, vector qx, vector qw) {
    int Q = rows(qx);
    vector[2] out;
    real ts = -log(n);
    real zl0 = (ts - 12 - t) / se;
    real zu0 = (ts + 6 - t) / se;
    if (zl0 >= 8) {
      // n p < exp(-12) wherever z < 8: first order, E[n p; z < 8]
      real a = n * exp(t + 0.5 * square(se) + std_normal_lcdf(8 - se | ));
      out[1] = -a; out[2] = -a;
      return out;
    }
    real zU = fmin(zu0, 8);
    real zL = zu0 <= -8 ? zl0 : fmax(zl0, -8);
    real h = (zU - zL) / 2;
    vector[Q] z = h * qx + (zU + zL) / 2;
    vector[Q] lw = log(h) + log(qw) - 0.5 * square(z) - 0.5 * log(2 * pi());
    vector[Q] lp = log_inv_logit(t + se * z);
    vector[Q] lf = n * log1m_inv_logit(t + se * z);
    real lJ = log_sum_exp(lw + lf);
    real lJd = log_sum_exp(lw + lf + lp);            // log E[f p] on the window
    if (zL == zl0 && zL > -30) {
      // closed-form left tail, first order in n p: E[1 - n p; z < zL] and E[p; z < zL]
      real lt = std_normal_lcdf(zL | );
      real lt2 = log(n) + t + 0.5 * square(se) + std_normal_lcdf(zL - se | );
      lJ = log_sum_exp(lJ, lt + log1m_exp(lt2 - lt));
      lJd = log_sum_exp(lJd, lt2 - log(n));
    }
    out[1] = lJ;
    out[2] = -n * exp(lJd - lJ);
    return out;
  }

  matrix logJ_table(vector tg, real n, real se, vector qx, vector qw) {
    int G = rows(tg);
    matrix[G, 2] LD;
    for (k in 1:G) LD[k] = logJ2(tg[k], n, se, qx, qw)';
    return LD;
  }

  // cubic Hermite interpolation of log J from the table; first-order closed
  // form left of the table, linear extrapolation right of it.
  real logJ_interp(real t, vector tg, real h, matrix LD, real n, real se) {
    int G = rows(tg);
    if (t <= tg[1]) return -n * exp(t + 0.5 * square(se));
    if (t >= tg[G]) return LD[G, 1] + LD[G, 2] * (t - tg[G]);
    int lo = 1;
    int hi = G;
    while (hi - lo > 1) {
      int mid = (lo + hi) %/% 2;
      if (t >= tg[mid]) lo = mid; else hi = mid;
    }
    real s = (t - tg[lo]) / h;
    real s2 = square(s);
    real s3 = s2 * s;
    return (2 * s3 - 3 * s2 + 1) * LD[lo, 1] + (s3 - 2 * s2 + s) * h * LD[lo, 2]
           + (-2 * s3 + 3 * s2) * LD[hi, 1] + (s3 - s2) * h * LD[hi, 2];
  }

  // log of the probability that a class-1 test with configuration linear
  // predictors eta (without its own effects) shows no failure, integrating
  // its effect x ~ N(0, su^2) by windowed Gauss-Legendre with a closed-form
  // left tail.
  real logI(row_vector eta, real n, real su, real se, vector glx, vector glw,
            vector tg, real h, matrix LD) {
    int M = rows(glx);
    int C = cols(eta);
    real lA = log(n) + 0.5 * square(se) + log_sum_exp(eta);
    real xs = -lA;
    real xL0 = xs - 12;
    real xU0 = xs + 6 + 4 * se;
    if (xL0 >= 8 * su) return 0;        // no failure is (almost) certain
    if (xU0 <= -8 * su) return -1e4;    // a failure is (almost) certain
    real xL = fmax(xL0, -8 * su);
    real xU = fmin(xU0, 8 * su);
    real hw = (xU - xL) / 2;
    vector[M] x = hw * glx + (xU + xL) / 2;
    vector[M] lf;
    for (m in 1:M) {
      real s = 0;
      for (c in 1:C) s += logJ_interp(eta[c] + x[m], tg, h, LD, n, se);
      lf[m] = s;
    }
    real lmid = log_sum_exp(log(hw) + log(glw) + lf
                            - 0.5 * square(x / su) - log(su) - 0.5 * log(2 * pi()));
    if (xL0 > -8 * su) {
      real lt = std_normal_lcdf(xL / su | );
      real lt2 = lA + 0.5 * square(su) + std_normal_lcdf((xL - square(su)) / su | );
      return log_sum_exp(lmid, lt + log1m_exp(lt2 - lt));
    }
    return lmid;
  }

  // log of int Binomial(y | n, inv_logit(eta + se z)) N(z | 0, 1) dz by
  // adaptive Gauss-Hermite quadrature: Newton steps for the mode on the
  // logit scale, then nodes centred there and scaled by the curvature.
  real log_binom_agq(int y, int n, real eta, real se, vector hx, vector hw) {
    int K = rows(hx);
    real s2 = square(se);
    real ph = (y + 0.5) / (n + 1.0);
    real ih = n * ph * (1 - ph);                   // data information at its mode
    real q = (ih * logit(ph) + eta / s2) / (ih + 1 / s2);   // normal-approximation start
    for (it in 1:10) {
      real p = inv_logit(q);
      real g = y - n * p - (q - eta) / s2;
      real h = -n * p * (1 - p) - 1 / s2;
      q -= fmax(fmin(g / h, 3), -3);    // damped Newton; the objective is concave
    }
    real p = inv_logit(q);
    real sq = 1 / sqrt(n * p * (1 - p) + 1 / s2);   // SD of the Laplace approximation, logit scale
    vector[K] qk = q + sq * hx;
    vector[K] lg = lchoose(n, y) + y * log_inv_logit(qk) + (n - y) * log1m_inv_logit(qk)
                   - 0.5 * square((qk - eta) / se) - log(se) - 0.5 * log(2 * pi())
                   + 0.5 * square(hx) + 0.5 * log(2 * pi());
    return log(sq) + log_sum_exp(log(hw) + lg);
  }

}
data {
  int<lower=0, upper=1> grouped;        // 0: one group, no group terms
  int<lower=1> P;                       // projects
  int<lower=1> C;                       // configurations
  int<lower=0> T;                       // tests with at least one failure
  array[T] int<lower=1, upper=P> proj;  // project of each such test
  array[P] int<lower=0, upper=1> cp;    // 1: centred group effects (data-rich group)
  int<lower=0> N1;                      // their cells with failures
  array[N1] int<lower=1, upper=T> t1;
  array[N1] int<lower=1, upper=C> c1;
  array[N1] int<lower=1> y1;
  array[N1] int<lower=1> n1;            // must equal nz (balanced runs)
  int<lower=0> N0;                      // their cells without failures
  array[N0] int<lower=1, upper=T> t0;
  array[N0] int<lower=1, upper=C> c0;
  array[N0] int<lower=1> n0;            // must equal nz (balanced runs)
  array[P] int<lower=0> K;              // tests with no failure, per project
  real<lower=1> nz;                     // runs per cell for those tests
  int<lower=2> M;                       // outer Gauss-Legendre nodes
  vector[M] glx;
  vector[M] glw;
  int<lower=2> Q;                       // inner Gauss-Legendre nodes
  vector[Q] qlx;
  vector[Q] qlw;
  int<lower=2> H;                       // adaptive Gauss-Hermite nodes
  vector[H] hx;
  vector[H] hw;
  int<lower=2> G;                       // grid for log J (spacing hg)
  vector[G] tg;
  real<lower=0> hg;
  // priors
  real gamma_loc; real<lower=0> gamma_scale;
  real alpha_loc; real<lower=0> alpha_scale;
  real<lower=0> sd_scale;               // student_t(3, 0, sd_scale) on every SD
}
parameters {
  real gamma;
  real alpha;
  vector<lower=0>[2] sd_ab;             // SDs of a_p (class) and b_p (intensity)
  cholesky_factor_corr[2] L_ab;
  matrix[2, P] ab_raw;                  // (a_p, b_p) if cp[p], else standardized
  real<lower=0> sd_u;
  vector[T] theta;                      // alpha + b_p + u_t for tests with failures
  real<lower=0> sd_v;
  vector[C] v;
  real<lower=0> sd_w;
  matrix[P, C] w_raw;                   // w_pc if cp[p], else standardized
  real<lower=0> sd_e;
}
transformed parameters {
  matrix[2, P] ab;
  matrix[P, C] w;
  {
    matrix[2, 2] Lc = diag_pre_multiply(sd_ab, L_ab);
    for (p in 1:P) {
      ab[:, p] = grouped == 0 ? rep_vector(0, 2) : (cp[p] ? ab_raw[:, p] : Lc * ab_raw[:, p]);
      w[p] = grouped == 0 ? rep_row_vector(0, C) : (cp[p] ? w_raw[p] : sd_w * w_raw[p]);
    }
  }
}
model {
  vector[P] a = ab[1]';
  vector[P] b = ab[2]';
  matrix[P, C] cfg = rep_matrix(v', P) + w;          // configuration terms per project
  matrix[P, C] eta0 = rep_matrix(alpha + b, C) + cfg;  // for tests without failures

  // tests with failures are flaky-capable
  for (t in 1:T) target += log_inv_logit(gamma + a[proj[t]]);
  matrix[G, 2] LD = logJ_table(tg, nz, sd_e, qlx, qlw);
  for (i in 1:N1) {
    real eta = theta[t1[i]] + cfg[proj[t1[i]], c1[i]];
    if (y1[i] == n1[i])   // every trial an event: J of the mirrored rate (n1 = nz)
      target += logJ_interp(-eta, tg, hg, LD, nz, sd_e);
    else
      target += log_binom_agq(y1[i], n1[i], eta, sd_e, hx, hw);
  }
  for (i in 1:N0)
    target += logJ_interp(theta[t0[i]] + cfg[proj[t0[i]], c0[i]], tg, hg, LD, nz, sd_e);

  // tests without failures, collapsed per project
  for (p in 1:P) {
    if (K[p] > 0) {
      real li = logI(eta0[p], nz, sd_u, sd_e, glx, glw, tg, hg, LD);
      target += K[p] * log_sum_exp(log1m_inv_logit(gamma + a[p]),
                                   log_inv_logit(gamma + a[p]) + li);
    }
  }

  // hierarchy
  theta ~ normal(alpha + b[proj], sd_u);
  v ~ normal(0, sd_v);
  {
    matrix[2, 2] Lc = diag_pre_multiply(sd_ab, L_ab);
    for (p in 1:P) {
      if (cp[p]) {
        ab_raw[:, p] ~ multi_normal_cholesky(rep_vector(0, 2), Lc);
        w_raw[p] ~ normal(0, sd_w);
      } else {
        ab_raw[:, p] ~ std_normal();
        w_raw[p] ~ std_normal();
      }
    }
  }

  // priors
  gamma ~ normal(gamma_loc, gamma_scale);
  alpha ~ normal(alpha_loc, alpha_scale);
  sd_ab ~ student_t(3, 0, sd_scale);
  sd_u ~ student_t(3, 0, sd_scale);
  sd_v ~ student_t(3, 0, sd_scale);
  sd_w ~ student_t(3, 0, sd_scale);
  sd_e ~ student_t(3, 0, sd_scale);
  L_ab ~ lkj_corr_cholesky(2);
}
generated quantities {
  real rho = L_ab[2, 1];
}
