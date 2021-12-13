// C++ port of PROPACK

#include "complex8.hpp"

// Globals
integer_t nopx = 0, nreorth = 0, nreorthu = 0, nreorthv = 0, ndot = 0, nitref = 0,
  nrestart = 0, nbsvd = 0, nlandim = 0, nsing = 0;
real_t tmvopx = 0.0, tgetu0 = 0.0, tupdmu = 0.0, tupdnu = 0.0, tintv = 0.0, tlanbpro = 0.0,
  treorthu = 0.0, treorthv = 0.0, telru = 0.0, telrv = 0.0, tbsvd = 0.0, tnorm2 = 0.0,
  tlansvd = 0.0, tritzvec = 0.0, trestart = 0.0, treorth = 0.0, tdot = 0.0;


// must include after frontmatter in this order!
#include "polyfill.hpp"
#include "utilities.hpp"
#include "csafescal.hpp"
#include "creorth.hpp"
#include "cgetu0.hpp"

void clanbpro(const integer_t& m,
	      const integer_t& n,
	      const integer_t& k0,
	      const integer_t& k,
	      const APROD_t aprod,
	      complex_t** U,
	      const integer_t& ldu,
	      complex_t** V,
	      const integer_t& ldv,
	      real_t** B,
	      const integer_t& ldb,
	      real_t& rnorm,
	      real_t soption[3],
	      integer_t ioption[2],
	      real_t* swork,
	      complex_t* cwork,
	      integer_t* iwork,
	      complex_t* cparm,
	      integer_t* iparm,
	      integer_t& ierr) {

  // %------------%
  // | Parameters |
  // %------------%
  constexpr integer_t one_int = 1;
  constexpr real_t one = 1.0;
  constexpr real_t zero = 0.0;
  constexpr real_t FUDGE = 1.01;
  constexpr real_t kappa = 0.717;
  constexpr char c_char = 'c';
  constexpr char n_char = 'n';

  // -------------------- Here begins executable code ---------------------
  auto t1 = std::chrono::steady_clock::now();
  std::chrono::time_point<std::chrono::steady_clock> t2;

  // %---------------------------------%
  // | Set machine dependent constants |
  // %---------------------------------%
  constexpr real_t eps = std::numeric_limits<real_t>::epsilon();
  constexpr real_t eps34 = std::pow(eps, 3.0/4.0);
  const real_t epsn = std::max(m, n)*eps;
  const real_t epsn2 = std::sqrt(std::max(m, n))*eps;

  // Parameters and declarations for goto nonsense
  real_t anormest;
  real_t amax, alpha, beta, mumax, numax, a1, b1, nrm;
  complex_t s;
  logical_t force_reorth;
  integer_t j0;
  constexpr integer_t imu = 0;
  const integer_t inu = imu + k;
  constexpr integer_t iidx = 0;
  constexpr integer_t is = 0;

  // %------------------------%
  // | Set default parameters |
  // %------------------------%
  const real_t delta = soption[0] < zero ? std::sqrt(eps/k) : soption[0];
  const real_t eta = soption[1] < zero ? eps34/std::sqrt(k) : soption[1];
  logical_t full_reorth = (delta <= eta) || (delta == zero);
  real_t anorm;
  if (soption[2] > zero) {
    anorm = soption[2];
  }
  else if (k0 > 0) {
    anorm = slapy2(&B[0][0], &B[0][1]);
    if (anorm <= zero) {
      ierr = -1;
      goto l9999;
    }
    else {
      anorm = zero;
    }
  }
  ierr = 0;

  // %---------------------%
  // | Get starting vector |
  // %---------------------%
  if (rnorm == zero) {
    constexpr integer_t ntry = 3; // NBPM: 3 not changed
    cgetu0(n_char, m, n, k0, ntry, &U[0][k0], rnorm, U,
	   ldu, aprod, cparm, iparm, ierr, ioption[0], anormest,
	   cwork);
    anorm = std::max(anorm, anormest);
  }

  // %------------------------------%
  // | Set pointers into work array |
  // %------------------------------%
  std::fill(swork, swork+2*k+2, 0);
  std::fill(iwork, iwork+2*k+1, 0);
  std::fill(cwork, cwork+std::max(m, n), 0);

  // %---------------------------%
  // | Prepare Lanczos iteration |
  // %---------------------------%
  if (k0 == 0) {
    amax = zero;
    alpha = zero;
    beta = rnorm;
    force_reorth = false;

    // %---------------------------------------------------%
    // | Compute ||A x|| / ||x|| for a random vector x     |
    // | to make it less likely that ||A|| is grossly      |
    // | underestimated at the beginning of the iteration. |
    // %---------------------------------------------------%
    if (n > m) {
      cgetu0(n_char, m, n, 0, 1, &cwork[is], nrm, U, ldu, aprod, cparm, // NBPM: 0, 1 unchanged
	     iparm, ierr, ioption[0], anormest, &cwork[is+m]);
    }
    else {
      cgetu0(c_char, m, n, 0, 1, &cwork[is], nrm, V, ldv, aprod, cparm, // NBPM: 0, 1 unchanged
	     iparm, ierr, ioption[0], anormest, &cwork[is+n]);
    }
    anorm = std::max(anorm, FUDGE*anormest);
    // std::cout << "anorm = " << anorm << std::endl;
    j0 = 0;
    if (beta != zero) {
      csafescal(&m, &beta, &U[0][0]);
    }
    // std::cout << "U0 = ";
    // for (auto i = 0; i < m; ++i) std::cout << U[i, 0] << " ";
    // std::cout << std::endl;
    swork[imu] = one;
    swork[inu] = one;
  }
  else {
    force_reorth = true;
    alpha = B[k0-1][0];
    beta = rnorm;
    if ((k0 < k) && (beta*delta < anorm*eps)) {
      full_reorth = true;
      ierr = k0;
    }

    iwork[iidx] = 0;
    iwork[iidx+1] = k0-1;
    iwork[iidx+2] = k0;

    csscal(&m, &rnorm, &U[0][k0], &one_int);
    t2 = std::chrono::steady_clock::now();
    creorth(m, k0, U, ldu, &U[0][k0], rnorm, &iwork[iidx], kappa,
	    &cwork[is], ioption[0]);
    treorthu += (std::chrono::steady_clock::now() - t2).count();
    csafescal(&m, &rnorm, &U[0][k0]);
    sset_mu(k0, &swork[imu], &iwork[iidx], epsn2);
    sset_mu(k0, &swork[inu], &iwork[iidx], epsn2);
    beta = rnorm;

    // %--------------------------------------%
    // | Estimate ||B||_2^2 as ||B^T * B||_1  |
    // %--------------------------------------%
    B[k0-1][1] = beta;
    amax = zero;
    for(integer_t j = 0; j < k0; ++j) {
      amax = std::max({amax, B[j][0], B[j][1]});
      if (j == 0) {
	anorm = std::max(anorm, FUDGE*alpha);
      }
      else if (j == 1) {
	a1 = B[0][1]/amax;
	a1 = FUDGE*amax*std::sqrt(std::pow(B[0][0]/amax, 2) + std::pow(a1, 2) + B[1][0]/amax*a1);
	anorm = std::max(anorm, a1);
      }
      else {
	a1 = B[j-1][0]/amax;
	b1 = B[j-1][1]/amax;
	a1 = FUDGE*amax*std::sqrt(std::pow(a1, 2) + std::pow(b1, 2) +
				  a1*B[j-2][1]/amax + B[j][0]/amax*b1);
	anorm = std::max(anorm, a1);
      }
    }
    j0 = k0;
  }
  numax = zero;
  mumax = zero;

  // %-------------------------------------------%
  // | Start Lanczos bidiagonalization iteration |
  // %-------------------------------------------%
  for (integer_t j = j0; j < k; ++j) {

    // %---------------------------------------------%
    // | alpha_{j} v_{j} = A'*u_{j} - beta_{j} v_{j} |
    // %---------------------------------------------%
    t2 = std::chrono::steady_clock::now();
    aprod(&c_char, m, n, &U[0][j], &V[0][j], cparm, iparm);
    tmvopx += (std::chrono::steady_clock::now() - t2).count();
    ++nopx;
    // std::cout << "V1 = ";
    // for (auto ii = 0; ii < m; ++ii) std::cout << V[ii, j] << " ";
    // std::cout << std::endl;

    if (j == 1) {
      alpha = scnrm2(&n, &V[0][j], &one_int);
      anorm = std::max(anorm, FUDGE*alpha);
      // std::cout << "j, alpha = " << j << ", " << alpha << std::endl;
    }
    else {
      const auto negbeta = -beta;
      csaxpy(&n, &negbeta, &V[0][j-1], &one_int, &V[0][j], &one_int);
      alpha = scnrm2(&n, &V[0][j], &one_int);
      // std::cout << "j, alpha = " << j << ", " << alpha << std::endl;

      // %------------------------------------%
      // | Extended local reorthogonalization |
      // %------------------------------------%
      t2 = std::chrono::steady_clock::now();
      if ((j > 0) && (ioption[1] > 0) && (alpha < kappa*beta)) {
	for (integer_t i = 0; i < ioption[1]; ++i) {
	  s = -cdotc(&n, &V[0][j-1], &one_int, &V[0][j], &one_int);
	  // std::cout << "s = " << s << std::endl;
	  caxpy(&n, &s, &V[0][j-1], &one_int, &V[0][j], &one_int);
	  nrm = scnrm2(&n, &V[0][j], &one_int);
	  // std::cout << "nrm = " << nrm << std::endl;
	  if (nrm >= kappa*alpha) {
	    goto l10;
	  }
	  alpha = nrm;
	}
      l10:
	swork[inu+j-2] = eps;
	alpha = nrm;
      }
      telrv += (std::chrono::steady_clock::now() - t2).count();

      B[j][0] = alpha;
      amax = std::max(amax, alpha);

      // %----------------------------%
      // | Update estimate of ||A||_2 |
      // %----------------------------%
      if (j == 1) {
	a1 = B[0][1]/amax;
	a1 = FUDGE*amax*std::sqrt(std::pow(B[0][0]/amax, 2) + std::pow(a1, 2) + B[1][0]/amax*a1);
      }
      else {
	a1 = B[j-1][0]/amax;
	b1 = B[j-1][1]/amax;
	a1 = FUDGE*amax*std::sqrt(std::pow(a1, 2) + std::pow(b1, 2) + a1*B[j-2][1]/amax + B[j][0]/amax*b1);
      }
      anorm = std::max(anorm, a1);
    }
    // std::cout << "j, alpha = " << j << ", " << alpha << std::endl;

    // %--------------------------%
    // | Update the nu recurrence |
    // %--------------------------%
    if (!full_reorth && (alpha != zero)) {
      supdate_nu(numax, &swork[imu], &swork[inu], j,
		 &B[0][0], &B[0][1], anorm, epsn2);
    }

    // %------------------------------%
    // | Reorthogonalize if necessary |
    // %------------------------------%
    if ((full_reorth || (numax > delta) || force_reorth) && (alpha != zero)) {
      //std::cout << "***REORTH***" << std::endl;
      if (full_reorth || (eta == zero)) {
	iwork[iidx] = 0;
	iwork[iidx+1] = j-1;
	iwork[iidx+2] = j;
      }
      else if (!force_reorth) {
	scompute_int(&swork[inu], j-1, delta, eta, &iwork[iidx]);
      }
      t2 = std::chrono::steady_clock::now();
      creorth(n, j-1, V, ldv, &V[0][j], alpha, &iwork[iidx],
	      kappa, &cwork[is], ioption[0]);
      treorthv += (std::chrono::steady_clock::now() - t2).count();

      sset_mu(j-1, &swork[inu], &iwork[iidx], eps);
      numax = eta;
      force_reorth = !force_reorth;
    }

    // %-----------------------------------------------%
    // | Check whether an invariant subspace was found |
    // %-----------------------------------------------%
    if ((alpha < anorm*epsn) && (j <= k)) {
      std::cerr << "UNVARIANT SUBSPACE, alpha small" << std::endl;
      rnorm = alpha;
      alpha = zero;
      // %------------------------------------------------%
      // | Try to build an orthogonal subspace, starting  |
      // | with a random vector.                          |
      // %------------------------------------------------%
      cgetu0(c_char, m, n, j-1, 3, &V[0][j], alpha, V, ldv, // NBPM: 3 left unchanged
	     aprod, cparm, iparm, ierr, ioption[0], anormest,
	     &cwork[is]);
      if (alpha == zero) {
	// %------------------------------------------------%
	// | We failed to generate a new random vector      |
	// | in span(A^H) orthogonal to span(V(:,1:j-1)).   |
	// | Most likely span(V(:,1:j-1)) is an invariant   |
	// | subspace.                                      |
	// %------------------------------------------------%
	// k = j; NBPM: this doesn't need to be changed
	ierr = -(j+1);
	goto l9999;
      }
      else {
	// %-------------------------------------------------%
	// | We have managed to generate a random vector     |
	// | in span(A^H) orthogonal to V(:,1:j-1), so we    |
	// | can continue the LBD and "deflate" the subspace |
	// | by setting alpha_{j} = 0.                       |
	// %-------------------------------------------------%
	csafescal(&n, &alpha, &V[0][j]);
	alpha = zero;
	force_reorth = (delta <= zero);
      }
    }
    else if ((j > 0) && !full_reorth && (j <= k) && (delta*alpha < anorm*eps)) {
      ierr = j+1;
    }
    B[j][0] = alpha;

    if (alpha != zero) {
      csafescal(&n, &alpha, &V[0][j]);
    }

    // std::cout << "V1_scaled = ";
    // for (auto ii = 0; ii < m; ++ii) std::cout << V[ii, j] << " ";
    // std::cout << std::endl;

    // %------------------------------------------------%
    // | beta_{j+1} u_{j+1} = A*v_{j} - alpha_{j} u_{j} |
    // %------------------------------------------------%
    t2 = std::chrono::steady_clock::now();
    aprod(&n_char, m, n, &V[0][j], &U[0][j+1], cparm, iparm);
    tmvopx += (std::chrono::steady_clock::now() - t2).count();
    ++nopx;

    const auto negalpha = -alpha;
    csaxpy(&m, &negalpha, &U[0][j], &one_int, &U[0][j+1], &one_int);
    beta = scnrm2(&m, &U[0][j+1], &one_int);
    // std::cout << "j, beta = " << j << ", " << beta << std::endl;

    // %------------------------------------%
    // | Extended local reorthogonalization |
    // %------------------------------------%
    t2 = std::chrono::steady_clock::now();
    if ((ioption[1] > 0) && (beta < kappa*alpha)) {
      for (integer_t i = 0; i < ioption[1]; ++i) {
	s = -cdotc(&m, &U[0][j], &one_int, &U[0][j+1], &one_int);
	// std::cout << "s = " << s << std::endl;
	caxpy(&m, &s, &U[0][j], &one_int, &U[0][j+1], &one_int);
	nrm = scnrm2(&m, &U[0][j+1], &one_int);
	// std::cout << "nrm = " << nrm << std::endl;
	if (nrm >= kappa*beta) {
	  goto l20;
	}
	beta = nrm;
      }
    l20:
      swork[imu+j-1] = eps;
      beta = nrm;
    }
    telru += (std::chrono::steady_clock::now() - t2).count();

    // std::cout << "j, beta = " << j << ", " << beta << std::endl;

    B[j][1] = beta;
    amax = std::max(amax, beta);

    // %----------------------------%
    // | Update estimate of ||A||_2 |
    // %----------------------------%
    if (j <= 0) {
      a1 = slapy2(&B[0][0], &B[0][1]);
    }
    else {
      a1 = B[j][0]/amax;
      a1 = amax*std::sqrt(std::pow(a1, 2) + std::pow(B[j][1]/amax, 2) + a1*B[j-1][1]/amax);
    }
    anorm = std::max(anorm, a1);

    // %--------------------------%
    // | Update the mu recurrence |
    // %--------------------------%
    if (!full_reorth && (beta != zero)) {
      supdate_mu(mumax, &swork[imu], &swork[inu], j, &B[0][0], &B[0][1], anorm, epsn2);
    }

    // %--------------------------------------%
    // | Reorthogonalize u_{j+1} if necessary |
    // %--------------------------------------%
    if ( (full_reorth || (mumax > delta) || force_reorth) && (beta != zero)) {
      if (full_reorth || (eta == zero)) {
	iwork[iidx] = 0;
	iwork[iidx+1] = j;
	iwork[iidx+2] = j+1;
      }
      else if (!force_reorth) {
	scompute_int(&swork[imu], j, delta, eta, &iwork[iidx]);
      }
      else {
	for (integer_t i = 0; i < 2*j+1; ++i) {
	  if (iwork[iidx+i-1] == j) {
	    iwork[iidx+i-1] = j+1;
	    goto l25;
	  }
	}
      }

    l25:
      t2 = std::chrono::steady_clock::now();
      creorth(m, j, U, ldu, &U[0][j+1], beta, &iwork[iidx], kappa, &cwork[is], ioption[0]);
      treorthu += (std::chrono::steady_clock::now() - t2).count();

      sset_mu(j, &swork[imu], &iwork[iidx], eps);
      mumax = eta;
      force_reorth = !force_reorth;
    }

    // %-----------------------------------------------%
    // | Check whether an invariant subspace was found |
    // %-----------------------------------------------%
    if ((beta < anorm*epsn) && (j <= k)) {
      std::cerr << "UNVARIANT SUBSPACE, beta small" << std::endl;
      rnorm = beta;
      beta = zero;
      // %-----------------------------------------------%
      // | Try to build an orthogonal subspace, starting |
      // | with a random vector.                         |
      // %-----------------------------------------------%
      cgetu0(n_char, m, n, j, 3, &U[0][j+1], beta, U, ldu, aprod, // NBPM: 3 left unchanged
	     cparm, iparm, ierr, ioption[0], anormest, &cwork[is]);
      if (beta == zero) {
	// %-----------------------------------------------%
	// | We failed to generate a new random vector     |
	// | in span(A) orthogonal to span(U(:,1:j)).      |
	// | Most likely span(U(:,1:j)) is an invariant    |
	// | subspace.                                     |
	// %-----------------------------------------------%
	// k = j+1; // NBPM: not necessary
	ierr = -(j+1);
	goto l9999;
      }
      else {
	// %------------------------------------------------%
	// | We have managed to generate a random vector    |
	// | in span(A) orthogonal to U(:,1:j), so we can   |
	// | continue the LBD and "deflate" the subspace by |
	// | setting beta_{j+1} = 0.                        |
	// %------------------------------------------------%
	csafescal(&m, &beta, &U[0][j+1]);
	beta = zero;
	force_reorth = (delta <= zero);
      }
    }
    else if (!full_reorth && (j <= k) && (delta*beta < anorm*eps)) {
      ierr = j+1;
    }

    B[j][1] = beta;
    if ((beta != zero) && (beta != one)) {
      csafescal(&m, &beta, &U[0][j+1]);
    }
    rnorm = beta;
    t2 = std::chrono::steady_clock::now();
  }
 l9999:
  soption[2] = anorm;
  tlanbpro += (std::chrono::steady_clock::now() - t1).count();
}


// **********************************************************************
