module KSmodel
  !ALL THE ROUTINES REQUIRED TO USE THE KS MODEL TURBULENT (LIKE) VELOCITY FIELD
  use cdata
  use general
    integer,private, parameter :: KSmodes=32
    real,private,dimension(3,KSmodes) :: unit_k,k,A,B
    real,private,dimension(3,Ksmodes) :: cross1,cross2
    real,private,dimension(KSmodes) :: omega
    contains 
subroutine get_KS_flow(x,u)
  
  ! Use our vectors to put together our KS velocity field
    implicit none
    real,intent(in)     :: x(3)
    real,intent(out)    :: u(3)
    real,dimension(3,KSmodes) :: cross1,cross2
    real :: argument
    real, dimension(3,KSmodes)::addition
    real :: sin_arg,cos_arg
    integer :: r
    
    u=0.
    argument=0.
    do r=1,KSmodes                 
    !find A(&B)crossed with unit_k as this appears in our sum for
    !the velocity field.
      cross1(1,r)=(A(2,r)*unit_k(3,r))-(A(3,r)*unit_k(2,r))
      cross1(2,r)=(A(3,r)*unit_k(1,r))-(A(1,r)*unit_k(3,r))
      cross1(3,r)=(A(1,r)*unit_k(2,r))-(A(2,r)*unit_k(1,r))
    
      cross2(1,r)=(B(2,r)*unit_k(3,r))-(B(3,r)*unit_k(2,r))
      cross2(2,r)=(B(3,r)*unit_k(1,r))-(B(1,r)*unit_k(3,r))
      cross2(3,r)=(B(1,r)*unit_k(2,r))-(B(2,r)*unit_k(1,r))
    
      argument=(k(1,r)*x(1))+(k(2,r)*x(2))+(k(3,r)*x(3))+(omega(r)*t)
   
      cos_arg=cos(argument)
      sin_arg=sin(argument)
      addition(1,r)=(cross1(1,r)*cos_arg)+(cross2(1,r)*sin_arg)
      addition(2,r)=(cross1(2,r)*cos_arg)+(cross2(2,r)*sin_arg)
      addition(3,r)=(cross1(3,r)*cos_arg)+(cross2(3,r)*sin_arg)
      u(1)=u(1)+addition(1,r)
      u(2)=u(2)+addition(2,r)
      u(3)=u(3)+addition(3,r)
    end do
  end subroutine get_KS_flow
    !************************************************************
    subroutine setup_KS
      implicit none

      real,dimension(3,KSmodes) :: orderK
      real,dimension(KSmodes) :: kk,delk,energy,klengths
      real :: arg
      real :: turn1,turnN
      integer :: i,j
    
    call wavenumbers(klengths)
  
    !put the wavelengths in order
    call order(klengths,KSmodes,1,kk) !the 1 means ascending order
  
    !sort the wave-vectors into the same order as their corresponding lengths
    do i=1,KSmodes
      do j=1,KSmodes
        if(kk(i)==klengths(j))then
          orderK(:,i)=k(:,j)
        end if
      end do
    end do
  
    k=orderK  !now we have k(3,N) that are listed in order of increasing length
  
    open(unit=19, file='./data/KSwavenumbers.log')
    do i=1,KSmodes
      unit_k(:,i)=k(:,i)/kk(i)
      write(19,*) kk(i), unit_k(:,i)
    end do
    close(19)    


      do i=1,KSmodes
      !now we find delk as defined in Malik & Vassilicos' paper
        if(i==1)delk(i)=(kk(i+1)-kk(i))*0.5
        if(i==KSmodes)delk(i)=(kk(i)-kk(i-1))*0.5
        if(i.gt.1.and.i.lt.KSmodes)delk(i)=(kk(i+1)-kk(i-1))*0.5
      end do
      
      !now find A&B that are perpendicular to each of our N wave-vectors
      call calc_KS_amplitudes(kk,delk,energy)
      
      do i=1,KSmodes
        arg=energy(i)*(kk(i)**3)            !these are used to define omega - the
        if(arg.gt.0.)omega(i)=sqrt(arg) !unsteadiness frequency (co-eff of t)
        if(arg==0.)omega(i)=0.
      end do
      
    ! Here we compute the turnover time of the fastest & slowest eddy
      turn1=2.0*pi/omega(1)   !the slowest one
      turnN=2.0*pi/omega(KSmodes)   !the fastest one
    
    ! Key information is output to a text file
      open(unit=57, file="./data/KSinfo.log")
        write(57,*) 'inner length scale',(2.*pi)/kk(KSmodes)
        write(57,*) 'outer length scale',(2.*pi)/kk(1)
        write(57,*) 'large eddie turnover time', turn1
        write(57,*) 'small eddie turnover time', turnN
        write(57,*) 'reynolds number',(kk(KSmodes)/kk(1))**(4./3.)  
      close(57)
      
    end subroutine
    !*****************************************************************
    subroutine wavenumbers(klengths)
    implicit none
    real,dimension(3) :: angle,dir_in
    integer, parameter :: total=100000
    real :: k_option(3,total), mkunit(total)
    real,dimension(KSmodes) :: klengths
    real:: bubble
    integer ::direction(3)
    integer :: i,j,num
    logical :: ne
    num=1
    !get the k-vectors now 
    do i=1,total  
      call random_number(angle)  
      !make sure none of the random numbers are zero at the given resolution 
      if(sum(angle)<epsilon(0.)) then
        call random_number(angle)
      end if
      angle=floor(9.*angle)/box_size
      call random_number(dir_in)
      direction=nint(dir_in)
      direction=2*direction -1  !positive or negative directions
       
      k_option(:,i)=direction(:)*pi*2.0*angle(:) !a possible wavevector
    
     !find the length of the current k_option vector
      mkunit(i)=sqrt((k_option(1,i)**2)+(k_option(2,i)**2)+(k_option(3,i)**2))
    
      if(i==1.and.mkunit(i).gt.0.)then ! first k_option vector is sure to be unique so we keep it
        k(:,num)=k_option(:,i)
        klengths(num)=mkunit(i)
      end if
    
      !now we check that the current length is unique (hasn't come before)
        if(i.gt.1.and.num.lt.KSmodes)then
          do j=i-1,1,-1
            if(mkunit(i).gt.0.0.and.mkunit(i) /= mkunit(j))then
                ne=.true.
            else
              ne=.false.
              exit
            end if
            if(j==1.and.ne)then !i.e. if length of current k_option is new...... 
              num=num+1
              k(:,num)=k_option(:,i) !load current k_option into k that we keep
              klengths(num)=mkunit(i)  ! store the length also
            end if
          end do
        end if     
      if(i==total.and.num.lt.KSmodes)print*,"Haven't got",KSmodes,"modes!!!!"
    end do
  end subroutine
    !****************************************************************************
    subroutine order(ad_, i_N, i_ord, B)
      !Bubble sorting algorithm
      implicit none
      integer, intent(in) :: i_N, i_ord
      real, intent(in)  :: ad_(i_N)
      real, intent(out) :: B(i_N)
      real :: c=1.
      integer :: n
      B = ad_
      do while(c>0.0)
         c = -1.0
         do n = 1, i_N-1
            if(   (i_ord==1 .and. B(n)>B(n+1))  &
             .or. (i_ord==2 .and. B(n)<B(n+1)) ) then
               c        = B(n)
               B(n)   = B(n+1)
               B(n+1) = c
               c = 1.
            end if
         end do
      end do
    end subroutine order
  !*************************************************************
  subroutine calc_KS_amplitudes(kk,delk,energy)
! Get our A's & B's that are perpendicular to k for each
! of our N wave-vectors k. This is done by making 2 extra random
! vectors, j & l and taking the cross product of these with k
    implicit none
    real,dimension(KSmodes),  intent(in)  :: kk,delk
    real,dimension(KSmodes),  intent(out) :: energy
    real :: j(3),l(3),newa(3),newa2(3)
    real :: ampA(KSmodes),ampB(KSmodes)
    integer :: r
    
    do r=1,KSmodes
       energy(r)=1.0+(kk(r))**2
       energy(r)=(kk(r)**2)*(energy(r)**(-11./6.)) !11./6. gives kolmogorov
       energy(r)=energy(r)*exp(-0.5*(kk(r)/kk(KSmodes))**2) 
       !set the lengths of A& B as defined in Malik & Vassilicos
       ampA(r)=sqrt(2.0*energy(r)*delk(r))!/3.0) 
       ampB(r)=ampA(r)                                
    
       call random_number(newa)
       call random_number(newa2)
       newa=2.*newa-1. ; newa2=2.*newa2-1.
       j=newa ; l=newa2  
       j=j/(sqrt(sum(newa**2)))
       l=l/(sqrt(sum(newa2**2)))
       
       !Now we call the subroutine that takes the vector product of k with
       !j (to get A) and with l (to get B)
       call get_vectorA(k,j,r,A)
       call get_vectorA(k,l,r,B)
     
      !Now that we have our unit A's & B's we multiply them by the amplitudes
      !we defined earlier.
      A(:,r)=ampA(r)*A(:,r)
      B(:,r)=ampB(r)*B(:,r)
      !These can be made to use cross_product from general.mod
      cross1(1,r)=(A(2,r)*unit_k(3,r))-(A(3,r)*unit_k(2,r))
      cross1(2,r)=(A(3,r)*unit_k(1,r))-(A(1,r)*unit_k(3,r))
      cross1(3,r)=(A(1,r)*unit_k(2,r))-(A(2,r)*unit_k(1,r))
  
      cross2(1,r)=(B(2,r)*unit_k(3,r))-(B(3,r)*unit_k(2,r))
      cross2(2,r)=(B(3,r)*unit_k(1,r))-(B(1,r)*unit_k(3,r))
      cross2(3,r)=(B(1,r)*unit_k(2,r))-(B(2,r)*unit_k(1,r))
    end do
  endsubroutine calc_KS_amplitudes
!***********************************************************************************
  subroutine get_vectorA(k,j,r,A)
  !aim:: to take the cross product of j(or l) with k to get a third
  !      vector that will be perpendicular to k. We then make this new
  !      vector have unit length so that we can multiply it with the
  !      correct amplitude later
    implicit none
    real,intent(in)  :: k(3,KSmodes),j(3)
    real,intent(out) :: A(3,KSmodes)
    real :: unity
    integer,intent(in) :: r
    !vector product (j X k)
    A(1,r)=(j(2)*k(3,r))-(j(3)*k(2,r))
    A(2,r)=(j(3)*k(1,r))-(j(1)*k(3,r))
    A(3,r)=(j(1)*k(2,r))-(j(2)*k(1,r))
    
    unity=sqrt((A(1,r)**2)+(A(2,r)**2)+(A(3,r)**2))
    
    !make A a unit vector
    A(:,r)=A(:,r)/unity
  end subroutine get_vectorA
end module
