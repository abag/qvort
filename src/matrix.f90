!>matrix all routines related to the matrix augmentation are in this module
module matrix
  contains
  !>get the inverse of n by n matrix
  subroutine findinv(matrix, inverse, n, errorflag)
	implicit none
	!declarations
	integer, intent(in) :: n
	integer, intent(out) :: errorflag  !return error status. -1 for error, 0 for normal
	real, intent(in), dimension(n,n) :: matrix  !input matrix
	real, intent(out), dimension(n,n) :: inverse !inverted matrix
	
	logical :: flag = .true.
	integer :: i, j, k, l
	real :: m
	real, dimension(n,2*n) :: augmatrix !augmented matrix
	
	!augment input matrix with an identity matrix
	do i = 1, n
		do j = 1, 2*n
			if (j <= n ) then
				augmatrix(i,j) = matrix(i,j)
			else if ((i+n) == j) then
				augmatrix(i,j) = 1
			else
				augmatrix(i,j) = 0
			endif
		end do
	end do
	
	!reduce augmented matrix to upper traingular form
	do k =1, n-1
		if (augmatrix(k,k) == 0) then
			flag = .false.
			do i = k+1, n
				if (augmatrix(i,k) /= 0) then
					do j = 1,2*n
						augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
					end do
					flag = .true.
					exit
				endif
				if (flag .eqv. .false.) then
					print*, "matrix is non - invertible"
					inverse = 0
					errorflag = -1
					return
				endif
			end do
		endif
		do j = k+1, n			
			m = augmatrix(j,k)/augmatrix(k,k)
			do i = k, 2*n
				augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
			end do
		end do
	end do
	
	!test for invertibility
	do i = 1, n
		if (augmatrix(i,i) == 0) then
			print*, "matrix is non - invertible"
			inverse = 0
			errorflag = -1
			return
		endif
	end do
	
	!make diagonal elements as 1
	do i = 1 , n
		m = augmatrix(i,i)
		do j = i , (2 * n)				
			   augmatrix(i,j) = (augmatrix(i,j) / m)
		end do
	end do
	
	!reduced right side half of augmented matrix to identity matrix
	do k = n-1, 1, -1
		do i =1, k
		m = augmatrix(i,k+1)
			do j = k, (2*n)
				augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
			end do
		end do
	end do				
	
	!store answer
	do i =1, n
		do j = 1, n
			inverse(i,j) = augmatrix(i,j+n)
		end do
	end do
	errorflag = 0
  end subroutine findinv
!subroutine to solve a set of simultaneous linear equations in 'n' variables
!using gaussian elimination method
!description:
!this subroutine finds the unknowns in a set of linear equations
!an example of using gaussian elimination is as follows:
!consider the equations:
!	      8x2 + 2x3 = -7
!	3x1 + 5x2 + 2x3 = 8
!	6x1 + 2x2 + 8x3 = 26
! write the equations in matrix form as 	0 8 2 -7
!                                               3 5 2 8
!                                               6 2 8 26
! re-order the equations as  6 2 8 26   
!                            3 5 2 8
!                            0 8 2 -7
! now convert to uppper traingular form to get   6 2 8 26
!                                                0 8 2 -7
!                                                0 0 -3 -3/2
! finally use back-substitution to solve the equation. for example, in case of x3 you get:
!                                                          -3x3 = -3/2
!                                                     =>     x3 = 1/2
!for a detailed explanation on the algorithm check out http://en.wikipedia.org/wiki gaussian_elimination                                                    
subroutine solve_linear_eqn(a, x, n, errflag)
!declarations
	implicit none
	integer, intent(in) :: n   !stores the number of unknowns
	integer, intent(out) :: errflag
	real, dimension(n,n+1) :: a!an n x n+1 matrix which stores the simultaneous equations
	real, intent(out), dimension(n) :: x !array to store the solutions
	integer :: i,j,k !counters for loops
	integer :: largest
	real :: temp, m, sums
	logical :: flag
	m = 0
	!solve using gaussian elimination
	flag = .false.
	x = 0
	do k = 1, n-1
		do j = 1, n
			if (a(j, 1) /= 0 ) flag = .true.
		end do
		if (flag .eqv. .false.) then
			print*,"no unique solution"
			errflag = -1
			x = 0
			exit
		else
			largest = k
			!find largest coefficient of first unknown
			do j = k, n
				if (abs(a(j, k)) > abs(a(largest,k))) largest = j
			end do			
			!make the equation with largest first coefficient as the first equation
			!largest coefficient is chosen to prevent round-off errors as far as possible
			do j = 1, n + 1
				temp = a(k, j)
				a(k,j) = a(largest,j)
				a(largest,j)=temp
			end do
		endif
		!convert the input matrix to upper traingualar form
		do j = k+1, n
			m = a(j,k)/a(k,k)
			do i = k+1, n+1
				a(j,i) = a(j,i) - m*a(k,i)
			end do
		end do
		!no unique solution exists if the last element in the upper triangular matrix is zero
		if (a(n,n) == 0) then
			print*,"no unique solution"
			errflag = -1
			x = 0
			exit
		else
			!find xn
			x(n) = a(n,n+1)/a(n,n)
			!find the remaining unknowns by back-substitution
			do i = n-1, 1, -1
				sums = 0
				do j = i+1, n
					sums = sums + a(i,j)*x(j)
				end do
				x(i) = (a(i,n+1) - sums)/a(i,i)
			end do
		endif
	end do
	
end subroutine
end module
