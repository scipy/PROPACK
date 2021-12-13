#ifndef UTILITIES_HPP_
#define UTILITIES_HPP_


int sgn(real_t val) {
    return (real_t(0) < val) - (val < real_t(0));
}


void sset_mu(const integer_t& k,
             real_t* mu,
             const integer_t* index,
             const real_t& val) {
  integer_t i = 0;
  integer_t p, q;
  while ((index[i] < k) && index[i] >= 0) {
    p = index[i];
    q = index[i+1];
    for (integer_t j = p; p < q; ++j) {
      mu[j] = val;
    }
    i += 2;
  }
}

// **********************************************************************

void scompute_int(const real_t* mu,
                  const integer_t& j,
                  const real_t& delta,
                  const real_t& eta,
                  integer_t* index) {

  auto t1 = std::chrono::steady_clock::now();
  if (delta < eta) {
    std::cerr << "Warning delta<eta in scompute_int" << std::endl;
    return;
  }

  integer_t s, k;
  integer_t ip = 0;
  index[0] = 0;
  integer_t i = 0;
  while(i <= j) {
    // find the next mu(k), k>i where abs(mu(k)) > delta
    for (k = i; k < j; ++k) {
      if (std::abs(mu[k]) > delta) {
        goto l10;
      }
    }
    goto l40;
    // find smallest i<k such that for all j=i,..,k, m(j) >= eta
  l10:
    for (s = k-1; s >= std::max(i, 0); --s) {
      if (std::abs(mu[s]) < eta) {
        goto l20;
      }
    }
  l20:
    ++ip;
    index[ip] = s + 1;
    for (i = s; i < j; ++i) {
      if (std::abs(mu[i]) < eta) {
        goto l30;
      }
    }
  l30:
    ++ip;
    index[ip] = i - 1;
  }
 l40:
  ++ip;
  index[ip] = j + 1;
  tintv += (std::chrono::steady_clock::now() - t1).count();
}

// **********************************************************************

void supdate_mu(real_t& mumax,
                real_t* mu,
                const real_t* nu,
                const integer_t& j,
                const real_t* alpha,
                const real_t* beta,
                const real_t& anorm,
                const real_t& eps1) {
  // %------------%
  // | Parameters |
  // %------------%
  constexpr real_t one = 1.0;

  // %--------------------%
  // | External Functions |
  // %--------------------%
  // real slamch
  // external slamch

  real_t d;
  auto t1 = std::chrono::steady_clock::now();
  if (j == 0) {
    d = eps1*(slapy2(&alpha[j], &beta[j]) + alpha[0]) + eps1*anorm;
    mu[0] = eps1/beta[0];
    mumax = std::abs(mu[0]);
  }
  else {
    mu[0] = alpha[0]*nu[0] - alpha[j]*mu[0];
    d = eps1*(slapy2(&alpha[j], &beta[j]) + alpha[0]) + eps1*anorm;
    mu[0] = (mu[0] + std::abs(d)*sgn(mu[0])) / beta[j];
    mumax = std::abs(mu[0]);
    for (integer_t k = 1; k < j-1; ++k) {
      mu[k] = alpha[k]*nu[k] + beta[k-1]*nu[k-1] - alpha[j]*mu[k];
      d = eps1*(slapy2(&alpha[j], &beta[j]) + slapy2(&alpha[k], &beta[k-1])) + eps1*anorm;
      mu[k] = (mu[k] + std::abs(d)*sgn(mu[k])) / beta[j];
      mumax = std::max(mumax, std::abs(mu[k]));
    }
    mu[j] = beta[j-1]*nu[j-1];
    d = eps1*(slapy2(&alpha[j], &beta[j]) + slapy2(&alpha[j], &beta[j-1])) + eps1*anorm;
    mu[j] = (mu[j] + std::abs(d)*sgn(mu[j])) / beta[j];
    mumax = std::max(mumax, std::abs(mu[j]));
  }
  mu[j+1] = one;
  tupdmu += (std::chrono::steady_clock::now() - t1).count();
}

// **********************************************************************

void supdate_nu(real_t& numax,
                const real_t* mu,
                real_t* nu,
                const integer_t j,
                const real_t* alpha,
                const real_t* beta,
                const real_t& anorm,
                const real_t& eps1) {

  // %------------%
  // | Parameters |
  // %------------%
  constexpr real_t one = 1.0;
  constexpr real_t zero = 0.0;

  real_t d;
  auto t1 = std::chrono::steady_clock::now();
  if (j > 0) {
    numax = zero;
    for (integer_t k = 0; k < j-1; ++k) {
      nu[k] = beta[k]*mu[k+1] + alpha[k]*mu[k] -beta[j-1]*nu[k];
      d = eps1*(slapy2(&alpha[k], &beta[k]) + slapy2(&alpha[j], &beta[j-1])) + eps1*anorm;
      nu[k] = (nu[k] + std::abs(d)*sgn(nu[k])) / alpha[j];
      numax = std::max(numax, std::abs(nu[k]));
    }
    nu[j] = one;
  }
  tupdnu += (std::chrono::steady_clock::now() - t1).count();
}

#endif // UTILITIES_HPP_
