!>all the routines required to setup and subsequently calculate the KS model
!>at a given position and time. This is used as as a normal fluid in the code.
!>The KS model is a model of a turbulent like velocity field created throught the
!>summation of random fourier modes
!!\f[
!!\mathbf{u}_n({\mathbf{x}},t)= \sum_{n=1}^{N_\textrm{KS}}\left(\mathbf{A}_n \times \mathbf{k}_n
!!\cos\phi_n + {\bf B}_n \times \mathbf{k}_n \sin\phi_n \right),
!!\f]
!!where \f$\phi_n=\mathbf{k}_n\cdot\mathbf{x}+\omega_n t\f$, \f$N_\textrm{KS}\f$ is the number of
!!modes, \f$\mathbf{k}_n\f$ and \f$\omega_n=k_n u_n\f$ are their wave vectors and frequencies.

module KSmodel
  use cdata
  use general
    !>arrays needed for the duration of the simualtion
    real,private,dimension(:,:),allocatable :: unit_k,k,A,B
    real,private,dimension(:,:),allocatable :: cross1,cross2
    real,private,dimension(:,:),allocatable:: addition !helper array for get_KS_flow
    real,private,dimension(:),allocatable :: omega !turnover time
    !>arrays needed for the setup only
    real,private,dimension(:,:),allocatable :: orderK
    real,private,dimension(:),allocatable  ::  klengths,kk,delk,KS_energy 
    contains 
    !************************************************************
    !>return the KS flow (u) at a given position x, t is a global
    !>variable which is used in here (flow is time dependent)
    subroutine get_KS_flow(x,u)
      implicit none
      real,intent(in)     :: x(3)
      real,intent(out)    :: u(3)
      real :: argument
      real :: sin_arg,cos_arg
      integer :: r
      u=0.
      argument=0.
      do r=1,KS_modes                 
        !find A(&B)crossed with unit_k as this appears in our sum for
        !the velocity field.
        argument=(k(1,r)*x(1))+(k(2,r)*x(2))+(k(3,r)*x(3))+(omega(r)*t)   
        cos_arg=cos(argument)
        sin_arg=sin(argument)

        addition(:,r)=(cross1(:,r)*cos_arg)+(cross2(:,r)*sin_arg)
        u(:)=u(:)+addition(:,r)
      end do
    end subroutine get_KS_flow
    !************************************************************
    !>the main routine to setup the KS flow, only performed once
    !>at the start of a run. KS vectors held within this module to
    !>be used by the get_KS_flow routine.
    subroutine setup_KS
      implicit none
      real :: arg
      real :: turn1,turnN
      integer :: i,j
      !begin by allocating the main KS arrays
      allocate(unit_k(3,KS_modes),k(3,KS_modes),A(3,KS_modes),B(3,KS_modes)) 
      allocate(cross1(3,KS_modes),cross2(3,KS_modes),addition(3,KS_modes),omega(KS_modes))
      !now allocate the arrays for this routine, they are deallocated at the end of the routine
      allocate(orderK(3,KS_modes))
      allocate(kk(KS_modes),delk(KS_modes),KS_energy(KS_modes),klengths(KS_modes))
      !if we have a mesh that we are printing to it is important the check the resolution
      if ((mesh_size>0).and.(mesh_size<4*KS_rey_int)) then
        call warning_message('ksmodel.mod, setup_KS','mesh resolution not sufficient to resolve KS model')
      end if
      call wavenumbers
  
      !put the wavelengths in order
      call order(klengths,KS_modes,1,kk) !the 1 means ascending order
  
      !sort the wave-vectors into the same order as their corresponding lengths
      do i=1,KS_modes
        do j=1,KS_modes
          if(kk(i)==klengths(j))then
            orderK(:,i)=k(:,j)
          end if
        end do
      end do
  
      k=orderK  !now we have k(3,N) that are listed in order of increasing length
  
      open(unit=19, file='./data/KSwavenumbers.log')
      do i=1,KS_modes
        unit_k(:,i)=k(:,i)/kk(i)
        write(19,*) kk(i), unit_k(:,i)
      end do
      close(19)    

      do i=1,KS_modes
      !now we find delk as defined in Malik & Vassilicos' paper
        if(i==1)delk(i)=(kk(i+1)-kk(i))*0.5
        if(i==KS_modes)delk(i)=(kk(i)-kk(i-1))*0.5
        if(i.gt.1.and.i.lt.KS_modes)delk(i)=(kk(i+1)-kk(i-1))*0.5
      end do      
      !now find A&B that are perpendicular to each of our N wave-vectors
      call calc_KS_amplitudes
      
      do i=1,KS_modes
        arg=KS_energy(i)*(kk(i)**3)            !these are used to define omega - the
        if(arg.gt.0.)omega(i)=sqrt(arg) !unsteadiness frequency (co-eff of t)
        if(arg==0.)omega(i)=0.
      end do
      
      !Here we compute the turnover time of the fastest & slowest eddy
      turn1=2.*pi/omega(1)   ; turnN=2.*pi/omega(KS_modes) 
    
      !Key information is output to a text file
      open(unit=57, file="./data/KSinfo.log")
        write(57,*) 'inner length scale',(2.*pi)/kk(KS_modes)
        write(57,*) 'outer length scale',(2.*pi)/kk(1)
        write(57,*) 'large eddie turnover time', turn1
        write(57,*) 'small eddie turnover time', turnN
        write(57,*) 'reynolds number',(kk(KS_modes)/kk(1))**(4./3.)  
      close(57)
      !also print to screen
      write(*,*) '-------------------KS info------------------'
      write(*,'(a,i5.4)') ' number of fourier modes:', KS_modes
      write(*,'(a,f6.2)') ' slope of spectrum:', KS_slope
      write(*,'(a,f6.2)') ' Reynolds number:', (kk(KS_modes)/kk(1))**(4./3.)
      write(*,'(a,f6.2,f6.2)') ' Large/small eddie turnover times:', turn1,turnN
      write(*,*) '--------------------------------------------'
      open(unit=77,file='./data/normal_timescale.log',status='replace')
        write(77,*) turn1 !write large eddy turnover as normal timescale
      close(77)
      !deallocate these helper arrays
      deallocate(orderK,kk,delk,KS_energy,klengths)
    end subroutine
    !*****************************************************************
    !>generate the KS wavenumbers
    subroutine wavenumbers
      implicit none
      real,dimension(3) :: angle,dir_in
      integer, parameter :: total=10000000
      real :: k_option(3,total), mkunit(total)
      real:: bubble
      integer ::direction(3)
      integer :: i,j,num
      logical :: ne
      !check the value of KS_rey_int to see if it will be displayed correctly
      num=1
      !get the k-vectors now 
      do i=1,total  
        call random_number(angle)  
        !make sure none of the random numbers are zero at the given resolution 
        if(sum(angle)<epsilon(0.)) then
          call random_number(angle)
        end if
        angle=floor((KS_rey_int+1)*angle)/box_size
        call random_number(dir_in)
        direction=nint(dir_in) ; direction=2*direction-1  !positive or negative directions
        if (KS_maximise_rey) then
          if (i==1) then
            write(*,*) 'minimising smallest wavenumber: may force anisotropy'
            angle(1)=0. ; angle(2)=0. ; angle(3)=1./box_size
            direction(3)=1
          end if
        end if
         
        k_option(:,i)=direction(:)*pi*2.0*angle(:) !a possible wavevector
      
        !find the length of the current k_option vector
        mkunit(i)=sqrt((k_option(1,i)**2)+(k_option(2,i)**2)+(k_option(3,i)**2))
        if(i==1.and.mkunit(i).gt.0.0)then
          k(:,num)=k_option(:,i)
          klengths(num)=mkunit(i)
        end if
      
        !now we check that the current length is unique (hasn't come before)
        if(i.gt.1.and.num.lt.KS_modes)then
          do j=i-1,1,-1
            if((mkunit(i)>0.).and.((mkunit(i)<mkunit(j)-0.001).or.(mkunit(i)>mkunit(j)+0.001)))then
                ne=.true.
            else
              ne=.false.
              exit
            end if
            if(j==1.and.ne)then !i.e. if length of current k_option is new
              num=num+1
              k(:,num)=k_option(:,i) !load current k_option into k that we keep
              klengths(num)=mkunit(i)  ! store the length also
            end if
          end do
        end if
        if (num==KS_modes) then
          write(*,*) ' found all KS wavenumbers'
          exit
        end if    
        if(i==total.and.num.lt.KS_modes) call fatal_error("ksmodel.mod","Haven't got enough modes to run KS")
      end do
    end subroutine
    !*************************************************************************
    !>bubble sorting algorithm used to order random KS modes
    subroutine order(ad_, i_N, i_ord, B)
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
    !>Get our A's & B's that are perpendicular to k for each
    !>of our N wave-vectors k. This is done by making 2 extra random
    !>vectors, j & l and taking the cross product of these with k
    subroutine calc_KS_amplitudes
      implicit none
      real :: j(3),l(3),newa(3),newa2(3)
      real,allocatable :: ampA(:),ampB(:)
      integer :: r 
      allocate(ampA(KS_modes),ampB(KS_modes))
      do r=1,KS_modes
         KS_energy(r)=1.0+(kk(r))**2
         KS_energy(r)=(kk(r)**2)*(KS_energy(r)**((KS_slope-2.)/2.)) !11./6. gives kolmogorov
         KS_energy(r)=KS_energy(r)*exp(-0.5*(kk(r)/kk(KS_modes))**2) 
         !set the lengths of A& B as defined in Malik & Vassilicos
         ampA(r)=sqrt(2.0*KS_energy(r)*delk(r))!/3.0) 
         ampB(r)=ampA(r)                                
      
         call random_number(newa) ; call random_number(newa2)
         newa=2.*newa-1. ; newa2=2.*newa2-1.
         j=newa ; l=newa2  
         j=j/(sqrt(sum(newa**2))) ; l=l/(sqrt(sum(newa2**2)))
         
         A(:,r)=cross_product(j(:),k(:,r))
         A(:,r)=A(:,r)/sqrt((A(1,r)**2)+(A(2,r)**2)+(A(3,r)**2)) 
         B(:,r)=cross_product(l(:),k(:,r))
         B(:,r)=B(:,r)/sqrt((B(1,r)**2)+(B(2,r)**2)+(B(3,r)**2))       
        !Now that we have our unit A's & B's we multiply them by the amplitudes
        !we defined earlier.
        A(:,r)=ampA(r)*A(:,r) ; B(:,r)=ampB(r)*B(:,r)
        !These can be made to use cross_product from general.mod
        cross1(:,r)=cross_product(A(:,r),unit_k(:,r))
        cross2(:,r)=cross_product(B(:,r),unit_k(:,r))
      end do
      deallocate(ampA,ampB)
    end subroutine
end module
