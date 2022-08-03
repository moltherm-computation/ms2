

module math_functions

  use ms2_global
  use math_types

  contains

    function addVectors(a, b) result(res)

      implicit none

      type(vector), intent(in) :: a, b
      type(vector) :: res

      res%x = a%x + b%x
      res%y = a%y + b%y
      res%z = a%z + b%z

    end function addVectors

    function subVectors(a, b) result(res)

      implicit none

      type(vector), intent(in) :: a, b
      type(vector) :: res

      res%x = a%x - b%x
      res%y = a%y - b%y
      res%z = a%z - b%z

    end function subVectors

    function scaleVector(a, factor) result(res)

      implicit none

      type(vector), intent(in) :: a
      real(RK), intent(in) :: factor
      type(vector) :: res

      res%x = factor * a%x
      res%y = factor * a%y
      res%z = factor * a%z

    end function scaleVector

end module math_functions
