#ifndef CREORTH_HPP_
#define CREORTH_HPP_

void ccgs(const integer_t& n,
          const integer_t& k,
          complex_t** V,
          const integer_t& ldv,
          complex_t* vnew,
          integer_t* index,
          complex_t* work) {

  // %-----------------%
  // | Local variables |
  // %-----------------%
  constexpr integer_t one_int = 1;
  constexpr complex_t zero_cmplx = complex_t(0.0, 0.0);
  constexpr complex_t one_cmplx = complex_t(1.0, 0.0);
  constexpr complex_t negone_cmplx = complex_t(-1.0, 0.0);
  integer_t i, p, q, l;
  auto ylocal = std::unique_ptr<complex_t[]>(new complex_t[n]);

  // The local variable ld was introduced to circumvent an apparent
  // bug in Intel's OpenMP implementation.
  const integer_t ld = ldv;
  constexpr integer_t tid = 0;
  constexpr integer_t nt = 1;
  integer_t cnk = n/nt;
  const integer_t st = tid*cnk + 1;

  i = 0;
  while((index[i] <= k) && (index[i] > 0)) {
    // Select the next block of columns from V
    p = index[i];
    q = index[i+1];
    l = q - p + 1;
    if (tid == 0) {
      ndot += l;
    }
    // Classical Gram-Schmidt: vnew = vnew - V(:,p:q)*(V(:,p:q)'*vnew)
    if (l > 0) {
      if (tid == nt - 1) {
        cnk = n-st+1;
      }
      cgemv('C', &cnk, &l, &one_cmplx, &V[st][p], &ld, &vnew[st], &one_int,
            &zero_cmplx, ylocal.get(), &one_int);

      if (tid == 0) {
        std::copy(ylocal.get(), ylocal.get() + n, work);
      }
      if (tid != 0) {
        std::transform(work, work + n, ylocal.get(), work, std::plus<>());
      }
      cgemv('N', &cnk, &l, &negone_cmplx, &V[st][p], &ld, work, &one_int,
            &zero_cmplx, ylocal.get(), &one_int);
      std::transform(&vnew[st-1], &vnew[st-1+cnk], ylocal.get(), &vnew[st-1], std::plus<>());
    }
    i += 2;
  }
}


void creorth(const integer_t& n,
             const integer_t& k,
             complex_t** V,
             const integer_t& ldv,
             complex_t* vnew,
             real_t& normvnew,
             integer_t* index,
             const real_t& alpha,
             complex_t* work,
             const integer_t& iflag) {

  // %------------%
  // | Parameters |
  // %------------%
  constexpr integer_t NTRY = 5;
  constexpr integer_t one_int = 1;

  // %-----------------%
  // | Local variables |
  // %-----------------%
  integer_t itry;
  real_t normvnew_0;

  if ((k <= 0) || (n <= 0)) {
    return;
  }
  auto t2 = std::chrono::steady_clock::now();

  for (integer_t itry = 0; itry < NTRY; ++itry) {
    normvnew_0 = normvnew;
    if (iflag == 1) {
      ccgs(n, k, V, ldv, vnew, index, work);
    } else {
      cmgs(n, k, V, ldv, vnew, index);
    }
    ndot += k;
    normvnew = scnrm2(&n, vnew, &one_int);
    if (normvnew > alpha*normvnew_0) {
      goto l9999_creorth;
    }
  }
  normvnew = 0.0;
  std::fill(vnew, vnew + n, 0);
 l9999_creorth:
  treorth += (std::chrono::steady_clock::now() - t2).count();
  nreorth += 1;
}

#endif // CREORTH_HPP_
