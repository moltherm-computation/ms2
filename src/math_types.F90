


module math_types

  use ms2_global

  type vector

    real(RK) :: x
    real(RK) :: y
    real(RK) :: z

  end type vector

  type matrix3D

    real(RK) :: entry(3, 3)

  end type matrix3D

end module math_types
