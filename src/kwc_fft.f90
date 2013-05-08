module kwc_fft
  use cdata
  use general
  implicit none
  integer*8 :: plan
  complex, allocatable :: in(:), out(:)
  complex, allocatable :: komega(:,:), komega_fft(:,:)
  real, allocatable :: k(:)
  integer :: komega_size, komega_shots, komega_count=0, komega_print=1
  contains
  subroutine initialise_fft
    implicit none
    real omega_max, omega_min, omega_track
    integer :: i
    allocate(in(pcount),out(pcount))
    allocate(k(pcount))
    do i=1,pcount
      if (i<=(pcount/2+1)) then
        k(i)=i-1
      else
        k(i)=-pcount/2+(i-pcount/2-1)
      end if
      !print*, k(i)
    end do
    k=k*2*pi/box_size
    !angular frequencies
    !omega_min=(log(1./(k(2)*corea)))*quant_circ*(k(2)**2)/(4*pi)
    !write(*,*), 'k min', k(2),'omega min', omega_min
    omega_min=(log(1./(k(2)*corea)))*quant_circ*(k(2)**2)/(4*pi)
    write(*,*), 'k min', k(2),'omega min', omega_min
    omega_track=(log(1./(k(5)*corea)))*quant_circ*(k(5)**2)/(4*pi)
    write(*,*), 'k track', k(5),'omega tracked', omega_track
    omega_max=(log(1./(k(pcount/2+1)*corea)))*quant_circ*(k(pcount/2+1)**2)/(4*pi)
    write(*,*)  'k max', k(pcount/2+1),'omega max', omega_max
    komega_shots=nint((5./omega_max)/dt)
    komega_size=nint((1./omega_track)/(5./omega_max))
    open(unit=33,file='./data/komega_dims.log')
      write(33,*), komega_shots, '%komega shots'
      write(33,*), komega_size, '%komega size'
    close(33)
    allocate(komega(komega_size,pcount))
  end subroutine
  subroutine perform_fft
   implicit none
   include 'fftw3.f'
   integer :: i
   character (len=40) :: print_file
   in=cmplx(f(:)%x(1),f(:)%x(2))

   !write(*,*) 'Input data for Fourier transform:'
   !do i=1,10
   !  write(*,*) in(i)
   !enddo

   ! Forward Fourier transform
   call dfftw_plan_dft_1d(plan,pcount,in,out,FFTW_FORWARD,FFTW_ESTIMATE)
   call dfftw_execute(plan)
   call dfftw_destroy_plan(plan)

   !Apply dissipation
   out(1)=0.
   out=out*exp(-(1E-8)*(k**4)*dt)
   out=out*exp(-(5)*(k**(-2))*dt)
   out(1)=0.
   !remove mean position
   if (mod(itime,shots)==0) then
      write(unit=print_file,fmt="(a,i4.4,a)")"./data/KWC_spect",itime/shots,".log"
      open(unit=37,file=print_file,status='replace')
       do i=1,pcount
         write(37,*) k(i), real(out(i)*conjg(out(i)))
       end do
     close(37)
   end if
     
   ! Inverse Fourier transform
   call dfftw_plan_dft_1d(plan,pcount,out,in,FFTW_BACKWARD,FFTW_ESTIMATE)
   call dfftw_execute(plan)
   call dfftw_destroy_plan(plan)
   in=in/pcount
   f(:)%x(1)=real(in)
   f(:)%x(2)=imag(in)
   if (mod(itime, komega_shots)==0) then
     komega_count=komega_count+1
     komega(komega_count,:)=in
     if (komega_count==komega_size) then
      write(unit=print_file,fmt="(a,i4.4,a)")"./data/komega",komega_print,".dat"
      open(unit=22,file=print_file,status='replace',form='unformatted',access='stream')
       write(22) komega 
      close(22)
      komega=0 ; komega_count=0 ; komega_print=komega_print+1
     end if
   end if
   !write(*,*) 'Recovered data from inverse Fourier transform:'
   !do i=1,10
   !  write(*,*) in(i)
   !enddo
  end subroutine
end module
