#ifndef CSAFESCAL_HPP_
#define CSAFESCAL_HPP_

void csafescal(const integer_t* n,
               const real_t* alpha,
               complex_t* x) {

  // %------------%
  // | Parameters |
  // %------------%
  constexpr integer_t one_int = 1;
  constexpr real_t one = 1.0;
  constexpr real_t zero = 0.0;

  // %-----------------%
  // | Local variables |
  // %-----------------%
  integer_t dummy, info;
  constexpr real_t sfmin = 1.0/std::numeric_limits<real_t>::max();

  // TODO: the right thing
  const real_t tmp = one / *alpha;
  csscal(n, &tmp, x, &one_int);

  // if (std::abs(alpha) >= sfmin) {
  //   const real_t tmp = one/alpha;
  //   csscal(n, &tmp, x, &one_int);
  // } else {
  //   clascl('General', &dummy, &dummy, alpha, &one, n, &one_int, x, n, &info);
  // }
}

#endif // CSAFESCAL_HPP
