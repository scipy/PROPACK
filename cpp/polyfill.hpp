#ifndef POLYFILL_HPP_
#define POLYFILL_HPP_

float slapy2(const float* x,
             const float* y) {
  return std::sqrt(*x**x + *y**y);
}

void csscal(const integer_t* n,
            const real_t* sa,
            complex_t* cx,
            /*unused*/const integer_t* incx) {
  struct scale {
    scale(const real_t* sa) : sa(*sa) { }
    complex_t operator()(complex_t c) { return sa*c; }
    real_t sa;
  };
  std::transform(cx, cx + *n, cx, scale(sa));
}


complex_t complex_prod(const complex_t& a, const complex_t& b) {
  return std::conj(a)*b;
}
complex_t add_c(const complex_t& a, const complex_t& b) {
  return a + b;
}
complex_t cdotc(const integer_t* n,
                const complex_t* cx,
                /*unused*/const integer_t* incx,
                const complex_t* cy,
                /*unused*/const integer_t* incy) {
  return std::inner_product(cx, cx + *n, cy, complex_t(), std::plus<>(), complex_prod);
}

template<typename VectorType, typename ScalarType>
void axpy(const integer_t* n,
          const ScalarType* a,
          const VectorType* x,
          /*unused*/const integer_t* incx,
          VectorType* y,
          /*unused*/const integer_t* incy) {
  // TODO: use std::transform
  for (std::size_t ii = 0; ii < *n; ++ii) {
    y[ii] += *a*x[ii];
  }
}

void caxpy(const integer_t* n,
           const complex_t* ca,
           const complex_t* cx,
           /*unused*/const integer_t* incx,
           complex_t* cy,
           /*unused*/const integer_t* incy) {
  return axpy(n, ca, cx, incx, cy, incy);
}

void csaxpy(const integer_t* n,
            const real_t* sa,
            const complex_t* cx,
            /*unused*/const integer_t* incx,
            complex_t* cy,
            /*unused*/const integer_t* incy) {
  return axpy(n, sa, cx, incx, cy, incy);
}

real_t scnrm2(const integer_t* n,
              const complex_t* cx,
              /*unused*/const integer_t* incx) {
  // TODO: use Boost norm functions
  real_t ans = 0;
  for (std::size_t ii = 0; ii < *n; ++ii) {
    ans += std::real(std::conj(cx[ii])*cx[ii]);
  }
  return std::sqrt(ans);
}

void clarnv(const integer_t* idist,
            /*unused*/integer_t* iseed,
            const integer_t* n,
            complex_t* x) {
  std::mt19937 gen;
  std::uniform_real_distribution<real_t> dis(-1.0, 1.0);

  switch(*idist) {
  case 1:
    // real and imaginary parts each uniform (0,1)
    break;
  case 2:
    // real and imaginary parts each uniform (-1,1)
    for (integer_t ii = 0; ii < *n; ++ii) {
      x[ii] = complex_t(dis(gen), dis(gen));
    }
    break;
  case 3:
    // real and imaginary parts each normal (0,1)
    break;
  case 4:
    // uniformly distributed on the disc abs(z) < 1
    break;
  case 5:
    // uniformly distributed on the circle abs(z) = 1
    break;
  }
}

#endif // POLYFILL_HPP_
