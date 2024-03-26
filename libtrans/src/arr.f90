program arr
   Implicit None
   Integer :: A(2, 3)
   A(:, :) = RESHAPE((/1, 2, 3, 4, 5, 6/), shape(A))
   write (*, '(2I)'), A(1, 3), A(2, 2)
end program arr
