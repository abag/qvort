module matrix
  contains
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
end module
