

!include this file with '#include "mathMacros.F90"'

!it is not possible to write e.g. SCALE_VECTOR(CROSS_PRODUCT(a, b), factor)

#define SQR(x) x*x

#define ADD_VECTOR(a, b) vector(a%x + b%x, a%y + b%y, a%z + b%z)
#define SUB_VECTOR(a, b) vector(a%x - b%x, a%y - b%y, a%z - b%z)
#define SCALE_VECTOR(a, factor) vector(factor * a%x, factor * a%y, factor * a%z)

#define SCALAR_PRODUCT(a, b) (a%x * b%x + a%y * b%y + a%z * b%z)

#define CROSS_PRODUCT(a, b) vector( \
    (a%y * b%z - a%z * b%y) , \
    (a%z * b%x - a%x * b%z) , \
    (a%x * b%y - a%y * b%x) )

#define SQR_VECTOR_NORM(a) SCALAR_PRODUCT(a, a)

#define SUB_ANINT_VECTOR(a) vector( \
    (a%x - anint(a%x)) , \
    (a%y - anint(a%y)) , \
    (a%z - anint(a%z)) )

#define VECTOR_MATRIX_PRODUCT(a, T) vector( \
    a%x * T%entry(1, 1) + a%y * T%entry(2, 1) + a%z * T%entry(3, 1), \
    a%x * T%entry(1, 2) + a%y * T%entry(2, 2) + a%z * T%entry(3, 2), \
    a%x * T%entry(1, 3) + a%y * T%entry(2, 3) + a%z * T%entry(3, 3) )







