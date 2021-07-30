#ifndef CGETU0_HPP_
#define CGETU0_HPP_

void cgetu0(const char transa,
            const integer_t& m,
            const integer_t& n,
            const integer_t& j,
            const integer_t& ntry,
            complex_t* u0,
            real_t& u0norm,
            complex_t** const U,
            const integer_t& ldu,
            APROD_t aprod,
            complex_t* dparm,
            integer_t* iparm,
            integer_t& ierr,
            const integer_t& icgs,
            real_t& anormest,
            complex_t* work) {

  // %------------%
  // | Parameters |
  // %------------%
  constexpr integer_t one_int = 1;
  constexpr real_t kappa = 0.717;

  // %-----------------%
  // | Local variables |
  // %-----------------%
  integer_t iseed[] = {1, 3, 5, 7};
  integer_t index[3];
  real_t nrm;

  // -------------------- Here begins executable code ---------------------
  auto t1 = std::chrono::steady_clock::now();
  std::chrono::time_point<std::chrono::steady_clock> t2;

  // %-------------------------%
  // | u0 is to be an m-vector |
  // | or an n-vector          |
  // %-------------------------%
  const bool is_mvec = (transa == 'n') || (transa == 'N');
  const integer_t rsize = is_mvec ? n : m;
  const integer_t usize = is_mvec ? m : n;

  const integer_t idist = 2; // real and imaginary parts each uniform (-1,1)
  ierr = 0;
  for (integer_t itry = 0; itry < ntry; ++itry) {
    clarnv(&idist, iseed, &rsize, work);
    nrm = scnrm2(&rsize, work, &one_int);
    t2 = std::chrono::steady_clock::now();
    aprod(&transa, m, n, work, u0, dparm, iparm);
    tmvopx += (std::chrono::steady_clock::now() - t2).count();
    ++nopx;

    u0norm = scnrm2(&usize, u0, &one_int);
    anormest = u0norm/nrm;

    if (j >= 0) {
      index[0] = 0;
      index[1] = j;
      index[2] = j+1;
      creorth(usize, j, U, ldu, u0, u0norm, index, kappa, work, icgs);
    }

    if (u0norm > 0) {
      goto l9999_cgetu0;
    }
  }
  ierr = -1;
 l9999_cgetu0:
  t2 = std::chrono::steady_clock::now();
  tgetu0 += (t2 - t1).count();
}

#endif // CGETU0_HPP_
