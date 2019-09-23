!> Module which contains initial conditions (lines) for filament (set in run.in) for more
!!information on the options see \ref INIT
module initial_line
  use cdata
  use general
  use periodic
  contains
  !****************************************************************************
  !>set up a single line from -z to z, number of particles is automatically adjusted
  !>to box size and \f$\delta\f$
  subroutine setup_single_line
    implicit none
    integer :: pcount_required
    integer :: mid_point
    integer :: i
    if (periodic_bc.or.periodic_bc_notx.or.periodic_bc_notxy) then
      !work out the number of particles required for single line
      !given the box size specified in run.i
      pcount_required=nint(box_size/(0.75*delta)) !75%
      write(*,*) 'changing size of pcount to fit with box_length and delta'
      write(*,*) 'pcount is now', pcount_required
      deallocate(f) ; pcount=pcount_required ; allocate(f(pcount))
    else
      call fatal_error('init.mod:setup_single_line', &
      'periodic boundary conditions required')
    end if
    mid_point=(pcount+1)/2
    do i=1, pcount
      f(i)%x(2)=0.!-box_size/3
      f(i)%x(1)=0.
      !f(i)%x(3)=-box_size/2.+box_size*real(i-1)/(pcount)
      f(i)%x(3)=-box_size/2.+box_size*real(2*i-1)/(2*pcount)
      if (i==1) then
        f(i)%behind=pcount ; f(i)%infront=i+1
      else if (i==pcount) then 
        f(i)%behind=i-1 ; f(i)%infront=1
      else
        f(i)%behind=i-1 ; f(i)%infront=i+1
      end if
      !zero the stored velocities
      f(i)%u1=0. ; f(i)%u2=0. ; f(i)%u3=0.
    end do
    if (kelvin_wave_casc) then
      !for kelvin wave cascade we allocate a new array
      allocate(KWC_z(pcount))
      KWC_z=f(:)%x(3)
    end if
  end subroutine
  !****************************************************************************
  !>set up a single line from -x to x, number of particles is automatically
  !adjusted
  !>to box size and \f$\delta\f$, small noise added to  y, z positions
  subroutine setup_xline_noise
    implicit none
    integer :: pcount_required
    integer :: i
    if (periodic_bc.or.periodic_bc_notx.or.periodic_bc_notxy) then
      !work out the number of particles required for single line
      !given the box size specified in run.i
      pcount_required=nint(box_size/(0.75*delta)) !75%
      write(*,*) 'changing size of pcount to fit with box_length and delta'
      write(*,*) 'pcount is now', pcount_required
      deallocate(f) ; pcount=pcount_required ; allocate(f(pcount))
    else
      call fatal_error('init.mod:setup_single_line', &
      'periodic boundary conditions required')
    end if
    do i=1, pcount
      f(i)%x(1)=-box_size/2.+box_size*real(2*i-1)/(2.*pcount)
      !add a little noise
      call random_number(f(i)%x(2:3))
      f(i)%x(2:3)=f(i)%x(2:3)*delta*0.5
      if (i==1) then
        f(i)%behind=pcount ; f(i)%infront=i+1
      else if (i==pcount) then
        f(i)%behind=i-1 ; f(i)%infront=1
      else
        f(i)%behind=i-1 ; f(i)%infront=i+1
      end if
      !zero the stored velocities
      f(i)%u1=0. ; f(i)%u2=0. ; f(i)%u3=0.
    end do
  end subroutine
  !****************************************************************************
  !>set up a single line from -z to z, number of particles is automatically adjusted
  !>to box size and \f$\delta\f$
  subroutine setup_soliton
    implicit none
    integer :: pcount_required
    integer :: i !for looping
    real :: ss !helper
    real :: theta, omega, soli_d !helpers
    if (periodic_bc.or.periodic_bc_notx.or.periodic_bc_notxy) then
      !work out the number of particles required for single line
      !given the box size specified in run.i
      pcount_required=nint(box_size/(0.6*delta)) !75%
      write(*,*) '2 soliton solution from Konno and Kakuhata (2005)'
      !--------------THESE ARE ALL SET IN RUN.IN---------------------
      write(*,'(a,f6.3,a,f6.3)') ' \lambda_r= ', soliton_lambda_r, ' \lambda_i= ', soliton_lambda_i
      write(*,'(a,f6.1,a,f6.1)') ' shift1= ', soliton_shift1, ' shift2= ', soliton_shift2
      write(*,'(a,f6.2,a,f6.2)') ' A1= ', soliton_A1, ' A2= ', soliton_A2
      !--------------------------------------------------------------
      write(*,*) 'changing size of pcount to fit with box_length and delta'
      write(*,*) 'pcount is now', pcount_required
      deallocate(f) ; pcount=pcount_required ; allocate(f(pcount))
    else
      call fatal_error('init.mod:setup_single_line', &
      'periodic boundary conditions required')
    end if
    soliton_lambda_i=soliton_lambda_i*2*pi/box_size !scale by box_size
    soliton_lambda_r=soliton_lambda_r*2*pi/box_size !scale by box_size
    soli_d=atan(-2*soliton_lambda_r*soliton_lambda_i/(soliton_lambda_i**2-soliton_lambda_r**2))
    do i=1, pcount
      ss=-box_size/2.+box_size*real(2*i-1)/(2.*pcount)
      f(i)%x(1)=soliton_A1*(soliton_lambda_i/(soliton_lambda_i**2+soliton_lambda_r**2))*&
                sin(2*soliton_lambda_r*(ss+soliton_shift1)+soli_d)*sech(2*soliton_lambda_i*(ss+soliton_shift1))+&
                soliton_A2*(soliton_lambda_i/(soliton_lambda_i**2+soliton_lambda_r**2))*&
                sin(2*soliton_lambda_r*(ss+soliton_shift2)+soli_d)*sech(2*soliton_lambda_i*(ss+soliton_shift2))
      f(i)%x(2)=-soliton_A1*(soliton_lambda_i/(soliton_lambda_i**2+soliton_lambda_r**2))*&
                cos(2*soliton_lambda_r*(ss+soliton_shift1)+soli_d)*sech(2*soliton_lambda_i*(ss+soliton_shift1))+&
                soliton_A2*(soliton_lambda_i/(soliton_lambda_i**2+soliton_lambda_r**2))*&
                cos(2*soliton_lambda_r*(ss+soliton_shift2)+soli_d)*sech(2*soliton_lambda_i*(ss+soliton_shift2))
      !f(i)%x(3)=-box_size/2.+box_size*real(2*i-1)/(2.*pcount)
      f(i)%x(3)=ss-soliton_A1*(soliton_lambda_i/(soliton_lambda_i**2+soliton_lambda_r**2))*&
                tanh(2*soliton_lambda_i*(ss+soliton_shift1))-&
                soliton_A2*(soliton_lambda_i/(soliton_lambda_i**2+soliton_lambda_r**2))*&
                tanh(2*soliton_lambda_i*(ss+soliton_shift2))
      if (i==1) then
        f(i)%behind=pcount ; f(i)%infront=i+1
      else if (i==pcount) then 
        f(i)%behind=i-1 ; f(i)%infront=1
      else
        f(i)%behind=i-1 ; f(i)%infront=i+1
      end if
      !zero the stored velocities
      f(i)%u1=0. ; f(i)%u2=0. ; f(i)%u3=0.
    end do
  end subroutine
  !****************************************************************************
  !>set up a helix from -z to z, number of particles is automatically adjusted
  !>to box size and \f$\delta\f$ 
  subroutine setup_helix
    implicit none
    integer :: pcount_required
    integer :: i
    if (periodic_bc.or.periodic_bc_notx.or.periodic_bc_notxy) then
      !work out the number of particles required for single line
      !given the box size specified in run.i
      pcount_required=nint(box_size/(0.65*delta)) !35%
      write(*,*) 'drawing a helix from -z to z'
      write(*,*) 'changing size of pcount to fit with box_length and delta'
      write(*,*) 'pcount is now', pcount_required
      deallocate(f) ; pcount=pcount_required ; allocate(f(pcount))
    else
      call fatal_error('init.mod:setup_single_line', &
      'periodic boundary conditions required')
    end if
    do i=1, pcount
      f(i)%x(3)=-box_size/2.+box_size*real(2*i-1)/(2.*pcount)
      f(i)%x(1)=5*delta*cos(2*pi*f(i)%x(3)/box_size)
      f(i)%x(2)=5*delta*sin(2*pi*f(i)%x(3)/box_size)
      if (i==1) then
        f(i)%behind=pcount ; f(i)%infront=i+1
      else if (i==pcount) then 
        f(i)%behind=i-1 ; f(i)%infront=1
      else
        f(i)%behind=i-1 ; f(i)%infront=i+1
      end if
      !zero the stored velocities
      f(i)%u1=0. ; f(i)%u2=0. ; f(i)%u3=0.
    end do
  end subroutine
!****************************************************************************
  !>hyperboloid - http://mathworld.wolfram.com/Hyperboloid.html
  subroutine setup_hyperboloid
    implicit none
    integer :: pcount_required
    integer :: line_size, line_position
    integer :: i, j
    real :: u, v !helper variables
    real :: shift1, shift2
    if (periodic_bc.or.periodic_bc_notx.or.periodic_bc_notxy) then
      !work out the number of particles required for single line
      !given the box size specified in run.in
      pcount_required=line_count*nint(box_size/(0.75*delta)) !75%
      print*, 'changing size of pcount to fit with box_length/delta and line_count'
      print*, 'pcount is now', pcount_required
      deallocate(f) ; pcount=pcount_required ; allocate(f(pcount))
    else
      call fatal_error('init.mod:setup_hyperboloid', &
      'periodic boundary conditions required')
    end if
    write(*,*) 'initf: hyperboloid,'
    write(*,'(a,f8.4)') ' width of bundle (at centre): ', hyperboloid_r*delta
    write(*,'(a,f8.4)') ' hyperboloid_e: ', hyperboloid_e
    line_size=int(pcount/line_count)
    do i=1, line_count
      v=pi*real(2*i-1)/line_count
      do j=1, line_size
        line_position=j+(i-1)*line_size
        u=-box_size/2.+box_size*real(2*j-1)/(2.*line_size)   
        f(line_position)%x(1)=hyperboloid_r*delta*sqrt(1.+(hyperboloid_e*u/box_size)**2)*cos(v) 
        f(line_position)%x(2)=hyperboloid_r*delta*sqrt(1.+(hyperboloid_e*u/box_size)**2)*sin(v)
        f(line_position)%x(3)=u
        if(j==1) then
          f(line_position)%behind=i*line_size
          f(line_position)%infront=line_position+1
        else if (j==line_size) then
          f(line_position)%behind=line_position-1
          f(line_position)%infront=(i-1)*line_size+1
        else
          f(line_position)%behind=line_position-1
          f(line_position)%infront=line_position+1
        end if
        f(line_position)%u1=0. ; f(line_position)%u2=0. ; f(line_position)%u3=0.
        if (j==1) then
          call random_number(shift1) ; call random_number(shift2)
          shift1=2*shift1-1 ; shift2=2*shift2-1
        end if
        f(line_position)%x(1)=f(line_position)%x(1)+2*shift1*delta
        f(line_position)%x(2)=f(line_position)%x(2)+2*shift2*delta
      end do
    end do
  end subroutine
  !*************************************************************************
  !>orthogonal lines from -z to z and -x to x which drives reconnection
  !>particle count automatically adjusted
  subroutine setup_orthog_lines
    implicit none
    integer :: pcount_required
    integer :: i 
    if (periodic_bc) then
      !work out the number of particles required for single line
      !given the box size specified in run.i
      pcount_required=2*nint(box_size/(0.75*delta)) !75%
      print*, 'changing size of pcount to fit with box_length and delta'
      print*, 'pcount is now', pcount_required
      deallocate(f) ; pcount=pcount_required ; allocate(f(pcount))
    else
      call fatal_error('init.mod:setup_crow', &
      'periodic boundary conditions required')
    end if
    write(*,*) 'initf: crow, separation of lines is:', 2.*delta 
    !loop over particles setting spatial and 'loop' position
    do i=1, pcount/2
      f(i)%x(1)=0.
      f(i)%x(3)=-box_size/2.+box_size*real(2*i-1)/(pcount)
      f(i)%x(2)=0.25*delta !-(delta/16.)*sin(pi*(box_size/2.+f(i)%x(3))/box_size)
      if (i==1) then
        f(i)%behind=pcount/2 ; f(i)%infront=i+1
      else if (i==pcount/2) then 
        f(i)%behind=i-1 ; f(i)%infront=1
      else
        f(i)%behind=i-1 ; f(i)%infront=i+1
      end if
      !zero the stored velocities
      f(i)%u1=0. ; f(i)%u2=0.; f(i)%u3=0.
    end do
    !second line
    do i=pcount/2+1, pcount
      f(i)%x(1)=box_size/2.-box_size*real(2*(i-pcount/2)-1)/(pcount)
      f(i)%x(2)=-0.25*delta!+(delta/16.)*sin(pi*(box_size/2.-f(i)%x(3))/box_size)
      f(i)%x(3)=0.
      if (i==(pcount/2+1)) then
        f(i)%behind=pcount ; f(i)%infront=i+1
      else if (i==pcount) then 
        f(i)%behind=i-1 ; f(i)%infront=1+pcount/2
      else
        f(i)%behind=i-1 ; f(i)%infront=i+1
      end if
      !zero the stored velocities
      f(i)%u1=0. ; f(i)%u2=0. ; f(i)%u3=0.
    end do    
  end subroutine
  !************************************************************************
  !>anti parallel lines from -z to z which drives the crow instability
  !>particle count automatically adjusted
  subroutine setup_crow
    implicit none
    integer :: pcount_required
    integer :: i 
    if (periodic_bc) then
      !work out the number of particles required for single line
      !given the box size specified in run.i
      pcount_required=2*nint(box_size/(0.75*delta)) !75%
      print*, 'changing size of pcount to fit with box_length and delta'
      print*, 'pcount is now', pcount_required
      deallocate(f) ; pcount=pcount_required ; allocate(f(pcount))
    else
      call fatal_error('init.mod:setup_crow', &
      'periodic boundary conditions required')
    end if
    write(*,*) 'initf: crow, separation of lines is:', 2.*delta 
    !loop over particles setting spatial and 'loop' position
    do i=1, pcount/2
      f(i)%x(1)=0.
      f(i)%x(3)=-box_size/2.+box_size*real(2*i-1)/(pcount)
      !tiny pertubation from straight line
      f(i)%x(2)=delta-(delta/16.)*sin(pi*(box_size/2.+f(i)%x(3))/box_size)
      !really bent lines below
      !f(i)%x(2)=16*delta-(12*delta)*sin(pi*(box_size/2.+f(i)%x(3))/box_size)
      !f(i)%x(2)=2*f(i)%x(2)
      if (i==1) then
        f(i)%behind=pcount/2 ; f(i)%infront=i+1
      else if (i==pcount/2) then 
        f(i)%behind=i-1 ; f(i)%infront=1
      else
        f(i)%behind=i-1 ; f(i)%infront=i+1
      end if
      !zero the stored velocities
      f(i)%u1=0. ; f(i)%u2=0.; f(i)%u3=0.
    end do
    !second line
    do i=pcount/2+1, pcount
      f(i)%x(1)=0.
      f(i)%x(3)=box_size/2.-box_size*real(2*(i-pcount/2)-1)/(pcount)
      !tiny pertubation from straight line
      f(i)%x(2)=-delta+(delta/16.)*sin(pi*(box_size/2.-f(i)%x(3))/box_size)
      !really bent lines below
      !f(i)%x(2)=-16*delta+(12*delta)*sin(pi*(box_size/2.-f(i)%x(3))/box_size)
      !f(i)%x(2)=2*f(i)%x(2)
      if (i==(pcount/2+1)) then
        f(i)%behind=pcount ; f(i)%infront=i+1
      else if (i==pcount) then 
        f(i)%behind=i-1 ; f(i)%infront=1+pcount/2
      else
        f(i)%behind=i-1 ; f(i)%infront=i+1
      end if
      !zero the stored velocities
      f(i)%u1=0. ; f(i)%u2=0. ; f(i)%u3=0.
    end do    
  end subroutine
  !*************************************************************************
  !> 4 bundles of lines in corners of the box
  subroutine setup_big_bundles
    implicit none
    integer :: pcount_required
    integer :: line_size, line_position
    integer :: i, j
    real :: rand1, rand2
    if (periodic_bc) then 
      !make sure line_count is a multiple of 4
      if (mod(line_count,4)/=0) then
        call fatal_error('init.mod:setup_big_bundles', &
        'line_count needs to be a multiple of 4')
       end if
      !work out the number of particles required for single line
      !given the box size specified in run.i
      pcount_required=line_count*nint(box_size/(0.75*delta)) !75%
      print*, 'changing size of pcount to fit with box_length/delta and line_count'
      print*, 'pcount is now', pcount_required
      deallocate(f) ; pcount=pcount_required ; allocate(f(pcount))
    else
      call fatal_error('init.mod:setup_big_bundles', &
      'periodic boundary conditions required')
    end if
    write(*,*) 'initf: big bundles' 
    line_size=int(pcount/line_count)
    do i=1, line_count
      rand1=runif(-1.,1.) ; rand2=runif(-1.,1.)
      do j=1, line_size
        line_position=j+(i-1)*line_size
        if (real(i)/line_count<=0.25) then
          f(line_position)%x(1)=-box_size/4.+(box_size/10.)*rand1
          f(line_position)%x(2)=-box_size/4.+(box_size/10.)*rand2
          f(line_position)%x(3)=-box_size/2.+box_size*real(2*j-1)/(2.*line_size)
        else if (real(i)/line_count<=0.5) then
          f(line_position)%x(1)=-box_size/4.+(box_size/10.)*rand1
          f(line_position)%x(2)=box_size/4.+(box_size/10.)*rand2
          f(line_position)%x(3)=box_size/2.-box_size*real(2*j-1)/(2.*line_size)
        else if (real(i)/line_count<=0.75) then
          f(line_position)%x(1)=box_size/4.+(box_size/10.)*rand1
          f(line_position)%x(2)=-box_size/4.+(box_size/10.)*rand2
          f(line_position)%x(3)=-box_size/2.+box_size*real(2*j-1)/(2.*line_size)
        else 
          f(line_position)%x(1)=box_size/4.+(box_size/10.)*rand1
          f(line_position)%x(2)=box_size/4.+(box_size/10.)*rand2
          f(line_position)%x(3)=box_size/2.-box_size*real(2*j-1)/(2.*line_size)
        end if
        if(j==1) then
          f(line_position)%behind=i*line_size
          f(line_position)%infront=line_position+1
        else if (j==line_size) then
          f(line_position)%behind=line_position-1
          f(line_position)%infront=(i-1)*line_size+1
        else
          f(line_position)%behind=line_position-1
          f(line_position)%infront=line_position+1
        end if
        f(line_position)%u1=0. ; f(line_position)%u2=0. ; f(line_position)%u3=0.
      end do
    end do
  end subroutine
!*************************************************************************
!> 4 bundles of lines in corners of the box
subroutine setup_bundles_recon
  implicit none
  integer :: pcount_required
  integer :: line_size, line_position
  integer :: i, j
  real :: rand1, rand2
  if (periodic_bc) then
    !make sure line_count is a multiple of 2
    if (mod(line_count,2)/=0) then
      call fatal_error('init.mod:setup_bundles_recon', &
      'line_count needs to be a multiple of 2')
     end if
    !work out the number of particles required for single line
    !given the box size specified in run.i
    pcount_required=line_count*nint(box_size/(0.75*delta)) !75%
    print*, 'changing size of pcount to fit with box_length/delta and line_count'
    print*, 'pcount is now', pcount_required
    deallocate(f) ; pcount=pcount_required ; allocate(f(pcount))
  else
    call fatal_error('init.mod:recon_bundles', &
    'periodic boundary conditions required')
  end if
  write(*,*) 'initf: recon bundles - bundle width', bundle_width, ' bundle separation', bundle_sep
  line_size=int(pcount/line_count)
  do i=1, line_count
    rand1=runif(0.,bundle_width) ; rand2=runif(0.,1.)
    do j=1, line_size
      line_position=j+(i-1)*line_size
      if (real(i)/line_count<=0.5) then
        f(line_position)%x(1)=bundle_sep/2+rand1*cos(rand2*2*pi)
        f(line_position)%x(2)=rand1*sin(rand2*2*pi)
        f(line_position)%x(3)=-box_size/2.+box_size*real(2*j-1)/(2.*line_size)
      else
        f(line_position)%x(1)=-bundle_sep/2+rand1*cos(rand2*2*pi)
        f(line_position)%x(2)=box_size/2.-box_size*real(2*j-1)/(2.*line_size)
        f(line_position)%x(3)=rand1*sin(rand2*2*pi)
      end if
      if(j==1) then
        f(line_position)%behind=i*line_size
        f(line_position)%infront=line_position+1
      else if (j==line_size) then
        f(line_position)%behind=line_position-1
        f(line_position)%infront=(i-1)*line_size+1
      else
        f(line_position)%behind=line_position-1
        f(line_position)%infront=line_position+1
      end if
      f(line_position)%u1=0. ; f(line_position)%u2=0. ; f(line_position)%u3=0.
    end do
  end do
  end subroutine
  !*************************************************************************
  !> 2 bundles of lines which reconnection
  subroutine setup_sultan_recon
    implicit none
    real, allocatable :: rad_shift(:), x_shift(:)
    real, allocatable :: rad_shift_r(:), rad_shift_theta(:)
    integer :: shift_help, shift_counter
    integer :: pcount_required, actual_line_count
    integer :: line_size, line_position
    integer :: i, j
    if (periodic_bc) then
      !work out the number of particles required for single line
      !given the box size specified in run.in
      actual_line_count=2*(3*line_count*(line_count-1)+1)
      pcount_required=actual_line_count*nint(box_size/(0.75*delta)) !75%
      print*, 'changing size of pcount to fit with box_length/delta and line_count'
      print*, 'pcount is now', pcount_required
      deallocate(f) ; pcount=pcount_required ; allocate(f(pcount))
    else
      call fatal_error('init.mod:setup_big_bundles', &
      'periodic boundary conditions required')
    end if
    write(*,*) 'initf: sultan reconnection'
    line_size=int(pcount/actual_line_count)
    allocate(rad_shift(actual_line_count/2),x_shift(actual_line_count/2))
    allocate(rad_shift_r(actual_line_count/2),rad_shift_theta(actual_line_count/2))
    rad_shift_r(1)=0. ; rad_shift_theta(1)=0.
    shift_counter=1
    if (line_count>1) then !we may not need this if statement
      do i=2,line_count
        do j=1,((i-1)*6)
          shift_counter=shift_counter+1
          rad_shift_r(shift_counter)=2*delta*(real(i)/line_count)
          rad_shift_theta(shift_counter)=2.*pi*real(j-1)/((i-1)*6)
        end do
      end do
    end if
    do i=1, actual_line_count/2
      do j=1, line_size
        line_position=j+(i-1)*line_size
        f(line_position)%x(1)=rad_shift_r(i)*cos(rad_shift_theta(i))-box_size/32
        f(line_position)%x(2)=rad_shift_r(i)*sin(rad_shift_theta(i))
        f(line_position)%x(3)=-box_size/2.+box_size*real(2*j-1)/(2.*line_size)
        if(j==1) then
          f(line_position)%behind=i*line_size
          f(line_position)%infront=line_position+1
        else if (j==line_size) then
          f(line_position)%behind=line_position-1
          f(line_position)%infront=(i-1)*line_size+1
        else
          f(line_position)%behind=line_position-1
          f(line_position)%infront=line_position+1
        end if
        f(line_position)%u1=0. ; f(line_position)%u2=0.
      end do
    end do
    do i=actual_line_count/2+1, actual_line_count
      do j=1, line_size
        line_position=j+(i-1)*line_size
        f(line_position)%x(1)=rad_shift_r(i-actual_line_count/2)*cos(rad_shift_theta(i-actual_line_count/2))+box_size/32
        f(line_position)%x(3)=rad_shift_r(i-actual_line_count/2)*sin(rad_shift_theta(i-actual_line_count/2))
        f(line_position)%x(2)=-box_size/2.+box_size*real(2*j-1)/(2.*line_size)
        if(j==1) then
          f(line_position)%behind=i*line_size
          f(line_position)%infront=line_position+1
        else if (j==line_size) then
          f(line_position)%behind=line_position-1
          f(line_position)%infront=(i-1)*line_size+1
        else
          f(line_position)%behind=line_position-1
          f(line_position)%infront=line_position+1
        end if
        f(line_position)%u1=0. ; f(line_position)%u2=0.
      end do
    end do
    deallocate(rad_shift,x_shift,rad_shift_r,rad_shift_theta)
  end subroutine
  !> 3 bundles of lines in all 3 cartesian directions
  subroutine setup_isotropic_bundles
    implicit none
    integer :: pcount_required
    integer :: line_size, line_position
    integer :: i, j
    real :: rand1, rand2
    if (periodic_bc) then 
      !make sure line_count is a multiple of 4
      if (mod(line_count,6)/=0) then
        call fatal_error('init.mod:setup_big_bundles', &
        'line_count needs to be a multiple of 6')
       end if
      !work out the number of particles required for single line
      !given the box size specified in run.i
      pcount_required=line_count*nint(box_size/(0.75*delta)) !75%
      print*, 'changing size of pcount to fit with box_length/delta and line_count'
      print*, 'pcount is now', pcount_required
      deallocate(f) ; pcount=pcount_required ; allocate(f(pcount))
    else
      call fatal_error('init.mod:setup_big_bundles', &
      'periodic boundary conditions required')
    end if
    write(*,*) 'initf: big bundles' 
    line_size=int(pcount/line_count)
    do i=1, line_count
      rand1=runif(-1.,1.) ; rand2=runif(-1.,1.)
      do j=1, line_size
        line_position=j+(i-1)*line_size
        !------------------lines -z to z--------------------
        if (real(i)/line_count<=0.166666666) then
          f(line_position)%x(1)=-box_size/4.+(box_size/10.)*rand1
          f(line_position)%x(2)=-box_size/4.+(box_size/10.)*rand2
          f(line_position)%x(3)=box_size/2.-box_size*real(2*j-1)/(2.*line_size)
        else if (real(i)/line_count<=0.3333333333) then
          f(line_position)%x(1)=box_size/4.+(box_size/10.)*rand1
          f(line_position)%x(2)=box_size/4.+(box_size/10.)*rand2
          f(line_position)%x(3)=-box_size/2.+box_size*real(2*j-1)/(2.*line_size)
        !------------------lines -y to y--------------------
        else if (real(i)/line_count<=0.5) then
          f(line_position)%x(1)=(box_size/10.)*rand1
          f(line_position)%x(2)=box_size/2.-box_size*real(2*j-1)/(2.*line_size)
          f(line_position)%x(3)=box_size/4.+(box_size/10.)*rand2
        else if (real(i)/line_count<=0.666666666) then
          f(line_position)%x(1)=(box_size/10.)*rand1
          f(line_position)%x(2)=-box_size/2.+box_size*real(2*j-1)/(2.*line_size)
          f(line_position)%x(3)=-box_size/4.+(box_size/10.)*rand2
        !------------------lines -x to x--------------------
        else if (real(i)/line_count<=0.833333333) then
          f(line_position)%x(1)=box_size/2.-box_size*real(2*j-1)/(2.*line_size)
          f(line_position)%x(2)=box_size/4.+(box_size/10.)*rand1
          f(line_position)%x(3)=box_size/4.+(box_size/10.)*rand2
        else 
          f(line_position)%x(1)=-box_size/2.+box_size*real(2*j-1)/(2.*line_size)
          f(line_position)%x(2)=-box_size/4.+(box_size/10.)*rand1
          f(line_position)%x(3)=-box_size/4.+(box_size/10.)*rand2
        end if
        if(j==1) then
          f(line_position)%behind=i*line_size
          f(line_position)%infront=line_position+1
        else if (j==line_size) then
          f(line_position)%behind=line_position-1
          f(line_position)%infront=(i-1)*line_size+1
        else
          f(line_position)%behind=line_position-1
          f(line_position)%infront=line_position+1
        end if
        f(line_position)%u1=0. ; f(line_position)%u2=0. ; f(line_position)%u3=0.
      end do
    end do
  end subroutine
  !*************************************************************************
  !> single bundle of lines in centre of box can be polarised or random
  !> depending on bundle_type parameter in run.in
  subroutine setup_central_bundle
    implicit none
    integer :: pcount_required
    integer :: line_size, line_position
    integer :: i, j
    real :: rand1, rand2, rand3
    real :: rr, theta !helper variables
    if (periodic_bc) then
      !make sure line_count is a multiple of 4
      !work out the number of particles required for single line
      !given the box size specified in run.i
      pcount_required=line_count*nint(box_size/(0.75*delta)) !75%
      print*, 'changing size of pcount to fit with box_length/delta and line_count'
      print*, 'pcount is now', pcount_required
      deallocate(f) ; pcount=pcount_required ; allocate(f(pcount))
    else
      call fatal_error('init.mod:setup_central_bundle', &
      'periodic boundary conditions required')
    end if
    write(*,*) 'initf: central bundle, bundle type: ', trim(bundle_type)
    select case(bundle_type)
      case('polarised')
        write(*,*) 'initf: central bundle, bundle distribution: ', trim(cent_bundle_dist)
    end select
    write(*,'(a,f8.4,a)') 'bundle will occupy ', lattice_ratio, ' of box' 
    line_size=int(pcount/line_count)
    do i=1, line_count
      rand1=runif(-1.,1.) ; rand2=runif(-1.,1.) !random bundle_type
      select case(cent_bundle_dist)
        case('uniform')
          theta=runif(0.,2*pi)
          rr=runif(0.,box_size*lattice_ratio)
        case('gaussian')
          theta=runif(0.,2*pi)
          rr=abs(rnorm(0.,(box_size*lattice_ratio/2.5)**2))
      end select
      do j=1, line_size
        line_position=j+(i-1)*line_size
        select case(bundle_type)
          case('polarised')
            f(line_position)%x(1)=rr*sin(theta)
            f(line_position)%x(2)=rr*cos(theta)
            f(line_position)%x(3)=-box_size/2.+box_size*real(2*j-1)/(2.*line_size)
          case('random')
            rand3=runif(0.,1.)
            if (rand3>0.5) then
              f(line_position)%x(1)=(box_size*lattice_ratio)*rand1/2.
              f(line_position)%x(2)=(box_size*lattice_ratio)*rand2/2.
              f(line_position)%x(3)=-box_size/2.+box_size*real(2*j-1)/(2.*line_size) 
            else
              f(line_position)%x(1)=(box_size*lattice_ratio)*rand1/2.
              f(line_position)%x(2)=-box_size/4.+(box_size*lattice_ratio)*rand2/2.
              f(line_position)%x(3)=box_size/2.-box_size*real(2*j-1)/(2.*line_size)
            end if
           case default
             call fatal_error('init.mod:setup_central_bundle', &
                  'bundle_type set to incorrect parameter')
        end select     
        if(j==1) then
          f(line_position)%behind=i*line_size
          f(line_position)%infront=line_position+1
        else if (j==line_size) then
          f(line_position)%behind=line_position-1
          f(line_position)%infront=(i-1)*line_size+1
        else
          f(line_position)%behind=line_position-1
          f(line_position)%infront=line_position+1
        end if
        f(line_position)%u1=0. ; f(line_position)%u2=0. ; f(line_position)%u3=0.
      end do
    end do
  end subroutine
  !*************************************************************************
  !>anti parallel lines from -z to z at -x, parallel lines from -z to z at x
  !>tests the smoothing routine, particle count automatically adjusted
  subroutine setup_smooth_test
    implicit none
    integer :: pcount_required
    integer :: i 
    if (periodic_bc) then
      !work out the number of particles required for single line
      !given the box size specified in run.i
      pcount_required=4*nint(box_size/(0.75*delta)) !75%
      print*, 'changing size of pcount to fit with box_length and delta'
      print*, 'pcount is now', pcount_required
      deallocate(f) ; pcount=pcount_required ; allocate(f(pcount))
    else
      call fatal_error('init.mod:setup_smooth_test', &
      'periodic boundary conditions required')
    end if
    write(*,*) 'initf: smooth_test, separation of lines is:', 2.*delta 
    !loop over particles setting spatial and 'loop' position   
    do i=1, pcount/4
      f(i)%x(1)=-box_size/4.
      f(i)%x(2)=delta
      f(i)%x(3)=-box_size/2.+2.*box_size*real(2*i-1)/(pcount)
      if (i==1) then
        f(i)%behind=pcount/4 ; f(i)%infront=i+1
      else if (i==pcount/4) then 
        f(i)%behind=i-1 ; f(i)%infront=1
      else
        f(i)%behind=i-1 ; f(i)%infront=i+1
      end if
      !zero the stored velocities
      f(i)%u1=0. ; f(i)%u2=0. ; f(i)%u3=0.
    end do
    !second line
    do i=pcount/4+1, 2*pcount/4
      f(i)%x(1)=-box_size/4.
      f(i)%x(2)=-delta
      f(i)%x(3)=-box_size/2.+2.*box_size*real(2*(i-pcount/4)-1)/(pcount)

      if (i==(pcount/4+1)) then
        f(i)%behind=2*pcount/4 ; f(i)%infront=i+1
      else if (i==2*pcount/4) then 
        f(i)%behind=i-1 ; f(i)%infront=1+pcount/4
      else
        f(i)%behind=i-1 ; f(i)%infront=i+1
      end if
      !zero the stored velocities
      f(i)%u1=0. ; f(i)%u2=0. ; ; f(i)%u3=0.
    end do    
    !third line
    do i=2*pcount/4+1, 3*pcount/4
      f(i)%x(1)=box_size/4.
      f(i)%x(2)=-delta
      f(i)%x(3)=-box_size/2.+2.*box_size*real(2*(i-pcount/2)-1)/(pcount)
      if (i==(2*pcount/4+1)) then
        f(i)%behind=3*pcount/4 ; f(i)%infront=i+1
      else if (i==3*pcount/4) then 
        f(i)%behind=i-1 ; f(i)%infront=1+2*pcount/4
      else
        f(i)%behind=i-1 ; f(i)%infront=i+1
      end if
      !zero the stored velocities
      f(i)%u1=0. ; f(i)%u2=0. ; f(i)%u3=0.
    end do
    !final line
    do i=3*pcount/4+1, pcount
      f(i)%x(1)=box_size/4.
      f(i)%x(2)=delta
      f(i)%x(3)=box_size/2.-2.*box_size*real(2*(i-3*pcount/4)-1)/(pcount)
      if (i==(3*pcount/4+1)) then
        f(i)%behind=pcount ; f(i)%infront=i+1
      else if (i==pcount) then 
        f(i)%behind=i-1 ; f(i)%infront=1+3*pcount/4
      else
        f(i)%behind=i-1 ; f(i)%infront=i+1
      end if
      !zero the stored velocities
      f(i)%u1=0. ; f(i)%u2=0. ; f(i)%u3=0.
    end do 
  end subroutine
  !*************************************************************************
  !>anti parallel lines from -z to z at -x, parallel lines from -z to z at x
  !>tests the smoothing routine, particle count automatically adjusted
  !>Kelvin waves added with imposed spectra
  subroutine setup_smooth_test_wave
    implicit none
    real :: wave_number, prefactor
    real :: amp, random_shift
    real :: xpos, ypos
    integer :: pcount_required
    integer :: line_size, line_position
    integer, parameter :: my_line_count=4
    integer :: i, j, k
    if (periodic_bc) then
      !work out the number of particles required for single line
      !given the box size specified in run.i
      pcount_required=my_line_count*nint(box_size/(0.55*delta)) !100% as waves are added
      write(*,*) 'changing size of pcount to fit with box_length and delta'
      write(*,'(a,i7.1)') ' pcount is now:', pcount_required
      deallocate(f) ; pcount=pcount_required ; allocate(f(pcount))
    else
      call fatal_error('init.mod:setup_smooth_test_wave', &
      'periodic boundary conditions required')
    end if
    write(*,'(a,i3.1,a)') ' drawing', my_line_count, ' lines from -z to +z'
    write(*,'(i4.1,a,a,a,f9.5)') wave_count, ' ',trim(wave_type),' wave pertubations, with spectral slope:', wave_slope
    line_size=int(pcount/my_line_count)
    !START THE LOOP
    do i=1, my_line_count
      if (i==1) then 
        xpos=-box_size/4.
        ypos=-2*delta
      else if (i==2) then
        xpos=-box_size/4.
        ypos=2*delta
      else if (i==3) then
        xpos=box_size/4.
        ypos=-2*delta
      else
        xpos=box_size/4.
        ypos=2*delta
      end if 
      do j=1, line_size
        line_position=j+(i-1)*line_size
        f(line_position)%x(1)=xpos
        f(line_position)%x(2)=ypos
        if (i==4) then
          f(line_position)%x(3)=-box_size/2.+box_size*real(2*j-1)/(2.*line_size)
        else
          f(line_position)%x(3)=box_size/2.-box_size*real(2*j-1)/(2.*line_size)        
        end if
        if(j==1) then
          f(line_position)%behind=i*line_size
          f(line_position)%infront=line_position+1
        else if (j==line_size) then
          f(line_position)%behind=line_position-1
          f(line_position)%infront=(i-1)*line_size+1
        else
          f(line_position)%behind=line_position-1
          f(line_position)%infront=line_position+1
        end if
        f(line_position)%u1=0. ; f(line_position)%u2=0. ; f(line_position)%u3=0.
      end do
      !we have now drawn the basic line, now we add the wave pertubations
      prefactor=wave_amp/(2.**wave_slope) !our starting wavenumber is 2
      if (i==1) then !only write this once
        write(*,*)'wave information recorded in ./data/wave_info.log'
        open(unit=34,file='./data/wave_info.log')
        write(34,*) '%------------------WAVE INFORMATION-------------------'
      end if
      do k=1, wave_count !wave_count set in run.in
        !wave_number=2+.1*k !starting wavenumber is 2, step in .1
        wave_number=3+.4*k !starting wavenumber is 2, step in .1
        amp=prefactor*(wave_number**wave_slope)
        call random_number(random_shift) !help things along with a 
        random_shift=random_shift*2*pi   !random shift \in (0,2\pi)
        if (i==1) then
          write(34,'(f9.5,f9.5,f9.5)') wave_number, amp, random_shift
        end if
        do j=1, line_size !reloop over particles
          line_position=j+(i-1)*line_size
          select case(wave_type) !what wave type do we want?
            case('planar')       !set this in run.in
              if (k==1) then !normal undefined for a straight line
                f(line_position)%x(1)=f(line_position)%x(1)+&
                amp*delta*sin(random_shift+wave_number*2.*pi*real(2.*j-1)/(2.*line_size))
              else
                f(line_position)%x(:)=f(line_position)%x(:)+&
                normalf(line_position)*amp*delta*sin(random_shift+wave_number*2.*pi*real(2.*j-1)/(2.*line_size))
              end if
            case('helical')
              if (k==1) then !normal/binormal undefined for a straight line
                f(line_position)%x(1)=f(line_position)%x(1)+&
                amp*delta*cos(random_shift+wave_number*2.*pi*real(2.*j-1)/(2.*line_size))
                f(line_position)%x(2)=f(line_position)%x(2)+&
                amp*delta*sin(random_shift+wave_number*2.*pi*real(2.*j-1)/(2.*line_size))
              else
                f(line_position)%x(:)=f(line_position)%x(:)+&
                normalf(line_position)*amp*delta*sin(random_shift+wave_number*2.*pi*real(2.*j-1)/(2.*line_size))+&
                binormalf(line_position)*amp*delta*cos(random_shift+wave_number*2.*pi*real(2.*j-1)/(2.*line_size))
              end if
            case default
              call fatal_error('initial.mod:smooth_wave_test','incorrect wave type parameter')
          end select
        end do
        call ghostp !we must call this routine at the end of adding every wave 
                    !the routines normalf, binormalf rely on correct ghostpoints
      end do !close k loop
      if (i==1) then
        close(34)
      end if 
    end do !closes the i loop
  end subroutine
  !*************************************************************************
  !>lines from the lop of the box to the bottom, helical/planar waves added with
  !>specific spectrum all set in run.in
  subroutine setup_wave_line
    implicit none
    real :: wave_number, prefactor
    real :: estimate_speed
    real :: amp, random_shift, random_shift2
    real :: xpos, ypos
    integer :: pcount_required
    integer :: line_size, line_position
    integer :: i, j, k
    !test run.in parameters, if wrong program will exit
    if (line_count==0) then
      call fatal_error('init.mod:setup_wave_line', &
      'you have not set a value for line_count in run.in')
    end if
    if (periodic_bc.or.periodic_bc_notx.or.periodic_bc_notxy) then
      !work out the number of particles required for single line
      !given the box size specified in run.i
      pcount_required=line_count*nint(box_size/(0.75*delta)) !100% as waves are added
      write(*,*) 'changing size of pcount to fit with box_length and delta'
      write(*,'(a,i7.1)') ' pcount is now:', pcount_required
      deallocate(f) ; pcount=pcount_required ; allocate(f(pcount))
    else
      call fatal_error('init.mod:setup_wave_line', &
      'periodic boundary conditions required')
    end if
    write(*,'(a,i3.1,a)') ' drawing', line_count, ' lines from -z to +z'
    write(*,'(i4.1,a,a,a,f9.5)') wave_count, ' ',trim(wave_type),' wave pertubations, with spectral slope:', wave_slope
    write(*,'(a,a)') ' wavenumber distribution:  ', wave_distribution
    write(*,'(a,i3.1)') ' starting wavenumber ', wave_start
    write(*,'(a,f9.5)') ' starting amplitude', wave_amp*delta
    write(*,'(a,i3.1)') ' wavenumber separation', wave_skip 
    if ((line_count==1).and.(wave_type=='helical')) then
      estimate_speed=abs((quant_circ*(2*pi*wave_start/box_size)**2/(4*pi))*(log(2./((2*pi*wave_start/box_size)*corea))-0.57721))
      !estimate_speed=(quant_circ/(2*pi*corea**2))*&
!(1.-sqrt(1+BesY0((2*pi*wave_start/box_size)*corea)/BesY1((2*pi*wave_start/box_size)*corea)))
      print*, 'just one helical wave, speed should be', estimate_speed
    end if
    line_size=int(pcount/line_count)
    !START THE LOOP
    do i=1, line_count
      if (line_count==1) then
        xpos=0. ; ypos=0.
      else
        call random_number(xpos) ; call random_number(ypos)
        xpos=(xpos-.5)*box_size*0.1 !1/10th of the box
        ypos=(ypos-.5)*box_size*0.1
      end if 
      do j=1, line_size
        line_position=j+(i-1)*line_size
        f(line_position)%x(1)=xpos
        f(line_position)%x(2)=ypos
        f(line_position)%x(3)=-box_size/2.+box_size*real(j-1)/line_size
        if(j==1) then
          f(line_position)%behind=i*line_size
          f(line_position)%infront=line_position+1
        else if (j==line_size) then
          f(line_position)%behind=line_position-1
          f(line_position)%infront=(i-1)*line_size+1
        else
          f(line_position)%behind=line_position-1
          f(line_position)%infront=line_position+1
        end if
        f(line_position)%u1=0. ; f(line_position)%u2=0. ; f(line_position)%u3=0.  
      end do
      !we have now drawn the basic line, now we add the wave pertubations
      prefactor=wave_amp/(wave_start**(wave_slope/2)) !our starting wavenumber is wave_start
      if (i==1) then !only write this once
        write(*,*)'wave information recorded in ./data/wave_info.log'
        open(unit=34,file='./data/wave_info.log')
        write(34,*) '%------------------WAVE INFORMATION-------------------'
      end if
      do k=1, wave_count !wave_count set in run.in
        select case(wave_distribution)
          case('spectrum')
            wave_number=wave_start+(k-1)*wave_skip
          case('gaussian')
             wave_number=nint(rnorm(real(wave_start),2.))
          case('bimodal')
             if (k<real(wave_count)/2) then
               wave_number=nint(rnorm(real(wave_start),0.))
             else
               wave_number=nint(rnorm(real(wave_start2),0.))
            end if
          case default
            call fatal_error('wave_line','incorrect parameter for wave_distribution')
        end select
        call random_number(random_shift) !help things along with a
        random_shift=random_shift*2*pi   !random shift \in (0,2\pi)
        call random_number(random_shift2) !help things along with a
        random_shift2=random_shift2*2*pi   !random shift \in (0,2\pi)
        if (wave_count==1) random_shift=0. !0 this for a single wave
        amp=prefactor*(wave_number**(wave_slope/2.))
        if (i==1) then
          write(34,'(f9.5,f9.5,f9.5)') wave_number, amp, random_shift
        end if
        do j=1, line_size !reloop over particles
          line_position=j+(i-1)*line_size
          select case(wave_type) !what wave type do we want?
            case('planar')       !set this in run.in
              if (k==1) then !normal undefined for a straight line
                f(line_position)%x(1)=f(line_position)%x(1)+&
                amp*delta*sin(random_shift+wave_number*2.*pi*real(2.*j-1)/(2.*line_size))
              else
                f(line_position)%x(:)=f(line_position)%x(:)+&
                normalf(line_position)*amp*delta*sin(random_shift+wave_number*2.*pi*real(2.*j-1)/(2.*line_size))
              end if
            case('helical')
              if (k==1) then !normal/binormal undefined for a straight line
                f(line_position)%x(1)=f(line_position)%x(1)+&
                amp*delta*cos(random_shift+wave_number*2.*pi*real(2.*j-1)/(2.*line_size))
                f(line_position)%x(2)=f(line_position)%x(2)+&
                amp*delta*sin(random_shift+wave_number*2.*pi*real(2.*j-1)/(2.*line_size))
              else
                f(line_position)%x(:)=f(line_position)%x(:)-&
                normalf(line_position)*amp*delta*cos(random_shift+wave_number*2.*pi*real(2.*j-1)/(2.*line_size))+&
                binormalf(line_position)*amp*delta*sin(random_shift+wave_number*2.*pi*real(2.*j-1)/(2.*line_size))
              end if
            case('KWC')
              f(line_position)%x(1)=f(line_position)%x(1)+&
              amp*delta*cos(random_shift+wave_number*2.*pi*real(j-1)/(line_size))
              f(line_position)%x(2)=f(line_position)%x(2)+&
              amp*delta*sin(random_shift+wave_number*2.*pi*real(j-1)/(line_size))
              f(line_position)%x(1)=f(line_position)%x(1)+&
              amp*delta*cos(random_shift2+wave_number*2.*pi*real(j-1)/(line_size))
              f(line_position)%x(2)=f(line_position)%x(2)-&
              amp*delta*sin(random_shift2+wave_number*2.*pi*real(j-1)/(line_size))
            case default
              call fatal_error('initial.mod:wave_line','incorrect wave type parameter')
          end select
        end do
        call ghostp !we must call this routine at the end of adding every wave 
                    !the routines normalf, binormalf rely on correct ghostpoints
      end do !close k loop
      if (i==1) then
        close(34)
      end if 
    end do !closes the i loop
    if (kelvin_wave_casc) then
      !for kelvin wave cascade we allocate a new array
      allocate(KWC_z(pcount))
      KWC_z=f(:)%x(3)
    end if
  end subroutine
  !*************************************************************************
  !>lines from the lop of the box to the bottom, helical/planar waves added with
  !>specific spectrum all set in run.in
  subroutine setup_hayder_wave
    implicit none
    real :: wave_number(5), amp(5), phase(5)
    integer :: pcount_required
    integer :: i, k
    !test run.in parameters, if wrong program will exit
    if (periodic_bc.or.periodic_bc_notx.or.periodic_bc_notxy) then
      !work out the number of particles required for single line
      !given the box size specified in run.i
      pcount_required=nint(box_size/(0.75*delta)) !100% as waves are added
      write(*,*) 'changing size of pcount to fit with box_length and delta'
      write(*,'(a,i7.1)') ' pcount is now:', pcount_required
      deallocate(f) ; pcount=pcount_required ; allocate(f(pcount))
    else
      call fatal_error('init.mod:setup_wave_line', &
      'periodic boundary conditions required')
    end if
    write(*,'(a)') ' drawing a line from -z to +z'
    write(*,'(a)') '5 helical wave pertubations '
    do i=1, pcount
      f(i)%x(1)=0. 
      f(i)%x(2)=0. 
      f(i)%x(3)=box_size/2.-box_size*real(2*i-1)/(2.*pcount)
      if (i==1) then
        f(i)%behind=pcount ; f(i)%infront=i+1
      else if (i==pcount) then 
        f(i)%behind=i-1 ; f(i)%infront=1
      else
        f(i)%behind=i-1 ; f(i)%infront=i+1
      end if
      !zero the stored velocities
      f(i)%u1=0. ; f(i)%u2=0. ; f(i)%u3=0.
    end do
    amp=(/5.5,4.5,3.8,4.0,4.5/)
    wave_number=(/5.,4.,3.,2.,1./)
    phase=(/.3422,.6343,.2213,.989342,.5345/)
    amp=amp*wave_amp
    phase=2*pi*phase
    do k=1, 5
      print*, k,' wavenumber ',wave_number(k),' amp ', amp(k), ' phase ', phase(k)
      do i=1, pcount
        f(i)%x(1)=f(i)%x(1)+&
                amp(k)*cos(wave_number(k)*2.*pi*real(2.*i-1)/(2.*pcount)+phase(k))
        f(i)%x(2)=f(i)%x(2)+&
                amp(k)*sin(wave_number(k)*2.*pi*real(2.*i-1)/(2.*pcount)+phase(k))
      end do 
    end do
  end subroutine
!*************************************************************************
  !>new initial condition
  subroutine setup_hayder_multi_wave
    implicit none
    real :: xpos, ypos
    real :: wave_number(5), amp(5), phase(line_count,5)
    integer :: pcount_required
    integer :: line_size, line_position
    integer :: i, j, k
    !test run.in parameters, if wrong program will exit
    if (line_count==0) then
      call fatal_error('init.mod:setup_wave_line', &
      'you have not set a value for line_count in run.in')
    end if
    if (periodic_bc.or.periodic_bc_notx.or.periodic_bc_notxy) then
      !work out the number of particles required for single line
      !given the box size specified in run.i
      pcount_required=line_count*nint(box_size/(0.75*delta)) !100% as waves are added
      write(*,*) 'changing size of pcount to fit with box_length and delta'
      write(*,'(a,i7.1)') ' pcount is now:', pcount_required
      deallocate(f) ; pcount=pcount_required ; allocate(f(pcount))
    else
      call fatal_error('init.mod:setup_wave_line', &
      'periodic boundary conditions required')
    end if
    write(*,'(a,i3.1,a)') ' drawing', line_count, ' lines from -z to +z'
    write(*,'(a)') ' adding 5 helical wave pertubations '
    line_size=int(pcount/line_count)
    !START THE LOOP
    do i=1, line_count
      if (i==1) then
        xpos=box_size/5. ; ypos=box_size/5.
      else if (i==2) then
        xpos=-box_size/5. ; ypos=-box_size/5.
      else
        call random_number(xpos) ; call random_number(ypos)
        xpos=(xpos-.5)*box_size*0.1 !1/10th of the box
        ypos=(ypos-.5)*box_size*0.1
      end if 
      do j=1, line_size
        line_position=j+(i-1)*line_size
        f(line_position)%x(1)=xpos
        f(line_position)%x(2)=ypos
        f(line_position)%x(3)=-box_size/2.+box_size*real(2*j-1)/(2.*line_size)
        if(j==1) then
          f(line_position)%behind=i*line_size
          f(line_position)%infront=line_position+1
        else if (j==line_size) then
          f(line_position)%behind=line_position-1
          f(line_position)%infront=(i-1)*line_size+1
        else
          f(line_position)%behind=line_position-1
          f(line_position)%infront=line_position+1
        end if
        f(line_position)%u1=0. ; f(line_position)%u2=0. ; f(line_position)%u3=0.  
      end do
      !---------wavenumbers---------------
      amp=(/5.5,4.5,3.8,4.0,4.5/)
      wave_number=(/5.,4.,3.,2.,1./)
      amp=amp*wave_amp
      !---------phases---------------
      call random_number(phase)
      phase=phase*2.*pi
      !phase(1,:)=0.
      do k=1, 5
        if (line_count==1) then
          print*, k,' wavenumber ',wave_number(k),' amp ', amp(k)
        end if
        do j=1, line_size !reloop over particles
          line_position=j+(i-1)*line_size
          f(line_position)%x(1)=f(line_position)%x(1)+&
          amp(k)*cos(phase(i,k)+wave_number(k)*2.*pi*real(2.*j-1)/(2.*line_size))
          f(line_position)%x(2)=f(line_position)%x(2)+&
          amp(k)*sin(phase(i,k)+wave_number(k)*2.*pi*real(2.*j-1)/(2.*line_size))
        end do
      end do !close k loop
    end do !closes the i loop
  end subroutine
  !*************************************************************************
  !>lines from the lop of the box to the bottom arranged in a lattice, the number of
  !>lines should be a square number
  subroutine setup_lattice
    implicit none
    real :: xpos, ypos
    integer :: pcount_required
    integer :: line_size, line_position
    integer :: counter=0
    integer :: i, j, k
    !test run.in parameters, if wrong program will exit
    if (line_count==0) then
      call fatal_error('init.mod:setup_lattice', &
      'you have not set a value for line_count in run.in')
    end if
    if (periodic_bc) then
      !work out the number of particles required for single line
      !given the box size specified in run.i
      pcount_required=line_count*nint(box_size/(.75*delta))
      write(*,*) 'changing size of pcount to fit with box_length and delta'
      write(*,'(a,i7.1)') ' pcount is now:', pcount_required
      deallocate(f) ; pcount=pcount_required ; allocate(f(pcount))
    else
      call fatal_error('init.mod:setup_lattice', &
      'periodic boundary conditions required')
    end if
    write(*,'(a,i3.1,a)') ' drawing', line_count, ' lines from -z to +z'
    write(*,'(a,f6.3,a)') ' lines in a lattice design occupying ', lattice_ratio, ' of box area'
    line_size=int(pcount/line_count)
    !START THE LOOP
    do i=1, floor(sqrt(real(line_count))) ; do k=1, floor(sqrt(real(line_count)))
      xpos=(-box_size/2.+box_size*((2.*i-1.)/(2*sqrt(real(line_count)))))*lattice_ratio
      ypos=(-box_size/2.+box_size*((2.*k-1.)/(2*sqrt(real(line_count)))))*lattice_ratio
      do j=1, line_size
        line_position=j+counter*line_size
        f(line_position)%x(1)=xpos
        f(line_position)%x(2)=ypos
        f(line_position)%x(3)=box_size/2.-box_size*real(2*j-1)/(2.*line_size)
        if(j==1) then
          f(line_position)%behind=(counter+1)*line_size
          f(line_position)%infront=line_position+1
        else if (j==line_size) then
          f(line_position)%behind=line_position-1
          f(line_position)%infront=counter*line_size+1
        else
          f(line_position)%behind=line_position-1
          f(line_position)%infront=line_position+1
        end if
        f(line_position)%u1=0. ; f(line_position)%u2=0. ; f(line_position)%u3=0.
      end do
      counter=counter+1
    end do ; end do
    if (counter/=line_count) then
      call fatal_error('init.mod:setup_lattice', &
      'line_count must be a square number')  
    end if
  end subroutine
  !*************************************************************************
  subroutine setup_tangle
    implicit none
    real :: rand1, rand2, rand3, rand4
    integer :: pcount_required
    integer :: line_size
    integer:: line_position
    integer :: i,j
    !test run.in parameters, if wrong program will exit
    if (line_count==0) then
      call fatal_error('init.mod:setup_tangle', &
      'you have not set a value for line_count in run.in')
    end if
    if (periodic_bc) then
      !work out the number of particles required for our lines
      !given the box size specified in run.i
      pcount_required=line_count*nint(box_size/(0.75*delta)) !75%
      write(*,*) 'changing size of pcount to fit with box_length and delta'
      write(*,*) 'pcount is now', pcount_required
      deallocate(f) ; pcount=pcount_required ; allocate(f(pcount))
    else
      call fatal_error('init.mod:setup_tangle',&
                       'periodic boundary conditions required for this initial condition') !cdata.mod
    end if
    line_size=int(pcount/line_count)
    do i=1, line_count
      call random_number(rand1)
      call random_number(rand2)
      call random_number(rand3)
      call random_number(rand4)
      rand1=(rand1-.5)*box_size
      rand3=(rand3-.5)*box_size
      if (rand4<0.5) then
        rand4=-1.
      else
        rand2=1.
      end if
      do j=1, line_size
        line_position=j+(i-1)*line_size
        if (rand2<0.5) then
          f(line_position)%x(1)=rand4*(box_size/10.)*sin(pi*real(2*j-1)/(2.*line_size))+rand3
          f(line_position)%x(2)=rand1
          f(line_position)%x(3)=-box_size/2.+box_size*real(2*j-1)/(2.*line_size)
        else
          f(line_position)%x(1)=rand1
          f(line_position)%x(2)=rand4*(box_size/10.)*sin(pi*real(2*j-1)/(2.*line_size))+rand3
          f(line_position)%x(3)=-box_size/2.+box_size*real(2*j-1)/(2.*line_size)
        end if
        if(j==1) then
          f(line_position)%behind=i*line_size
          f(line_position)%infront=line_position+1
        else if (j==line_size) then
          f(line_position)%behind=line_position-1
          f(line_position)%infront=(i-1)*line_size+1
        else
          f(line_position)%behind=line_position-1
          f(line_position)%infront=line_position+1
        end if
        f(line_position)%u1=0. ; f(line_position)%u2=0. ; f(line_position)%u3=0.
      end do
    end do
  end subroutine
  !*************************************************************************
  subroutine setup_line_of_lines
    implicit none
    real :: rand1, rand2, rand3, rand4
    real :: rand5, rand6
    integer :: pcount_required
    integer :: line_size
    integer:: line_position
    integer :: i,j
    !test run.in parameters, if wrong program will exit
    if (line_count==0) then
      call fatal_error('init.mod:setup_tangle', &
      'you have not set a value for line_count in run.in')
    end if
    if (periodic_bc) then
      !work out the number of particles required for our lines
      !given the box size specified in run.i
      pcount_required=line_count*nint(box_size/(0.75*delta)) !75%
      write(*,*) 'changing size of pcount to fit with box_length and delta'
      write(*,*) 'pcount is now', pcount_required
      deallocate(f) ; pcount=pcount_required ; allocate(f(pcount))
    else
      call fatal_error('init.mod:setup_criss_cross',&
                       'periodic boundary conditions required for this initial condition') !cdata.mod
    end if
    line_size=int(pcount/line_count)
    do i=1, line_count
        if (mod(line_count,2)==0) then
          rand2=1.
        else
          rand2=-1.
        end if
        rand1=-box_size/2+(real(i)/(line_count+1))*box_size
        do j=1, line_size
          line_position=j+(i-1)*line_size
          f(line_position)%x(1)=rand1
          f(line_position)%x(2)=0.
          f(line_position)%x(3)=rand2*(-box_size/2.+box_size*real(2*j-1)/(2.*line_size))
          if(j==1) then
            f(line_position)%behind=i*line_size
            f(line_position)%infront=line_position+1
          else if (j==line_size) then
            f(line_position)%behind=line_position-1
            f(line_position)%infront=(i-1)*line_size+1
          else
            f(line_position)%behind=line_position-1
            f(line_position)%infront=line_position+1
          end if
          f(line_position)%u1=0. ; f(line_position)%u2=0. ; f(line_position)%u3=0.
        end do
    end do
  end subroutine
  !*************************************************************************
  subroutine setup_criss_cross
    implicit none
    real :: rand1, rand2, rand3, rand4
    real :: rand5, rand6
    integer :: pcount_required
    integer :: line_size
    integer:: line_position
    integer :: i,j
    !test run.in parameters, if wrong program will exit
    if (line_count==0) then
      call fatal_error('init.mod:setup_tangle', &
      'you have not set a value for line_count in run.in')
    end if
    if (periodic_bc) then
      !work out the number of particles required for our lines
      !given the box size specified in run.i
      pcount_required=line_count*nint(box_size/(0.75*delta)) !75%
      write(*,*) 'changing size of pcount to fit with box_length and delta'
      write(*,*) 'pcount is now', pcount_required
      if (criss_cross_bundle>1) then
        write(*,'(a,i4.3)') ' putting in some bundling with size: ', criss_cross_bundle
        write(*,'(a,f6.4)') ' width of bundles ', criss_cross_width*delta
      end if
      deallocate(f) ; pcount=pcount_required ; allocate(f(pcount))
    else
      call fatal_error('init.mod:setup_criss_cross',&
                       'periodic boundary conditions required for this initial condition') !cdata.mod
    end if
    line_size=int(pcount/line_count)
    do i=1, line_count
      if ((mod(i,criss_cross_bundle)==0).or.(i==1)) then
        call random_number(rand1)
        call random_number(rand2)
        call random_number(rand3)
        call random_number(rand4)
        rand1=(rand1-.5)*box_size
        rand2=(rand2-.5)*box_size
        if (rand4<0.5) then
          rand4=-1
         else
          rand4=1
         end if
      end if
      rand5=runif(-criss_cross_width*delta,criss_cross_width*delta)
      rand6=runif(-criss_cross_width*delta,criss_cross_width*delta)
      if (rand3<0.33333333) then
        do j=1, line_size
          line_position=j+(i-1)*line_size
          f(line_position)%x(1)=rand1+rand5
          f(line_position)%x(2)=rand2+rand6
          f(line_position)%x(3)=rand4*(-box_size/2.+box_size*real(2*j-1)/(2.*line_size))
          if(j==1) then
            f(line_position)%behind=i*line_size
            f(line_position)%infront=line_position+1
          else if (j==line_size) then
            f(line_position)%behind=line_position-1
            f(line_position)%infront=(i-1)*line_size+1
          else
            f(line_position)%behind=line_position-1
            f(line_position)%infront=line_position+1
          end if
          f(line_position)%u1=0. ; f(line_position)%u2=0. ; f(line_position)%u3=0.
        end do
      else if (rand3<0.6666666) then
        do j=1, line_size
          line_position=j+(i-1)*line_size
          f(line_position)%x(1)=rand1+rand5
          f(line_position)%x(3)=rand2+rand6
          f(line_position)%x(2)=rand4*(-box_size/2.+box_size*real(2*j-1)/(2.*line_size))
          if(j==1) then
            f(line_position)%behind=i*line_size
            f(line_position)%infront=line_position+1
          else if (j==line_size) then
            f(line_position)%behind=line_position-1
            f(line_position)%infront=(i-1)*line_size+1
          else
            f(line_position)%behind=line_position-1
            f(line_position)%infront=line_position+1
          end if
          f(line_position)%u1=0. ; f(line_position)%u2=0. ; f(line_position)%u3=0.
        end do
      else
        do j=1, line_size
          line_position=j+(i-1)*line_size
          f(line_position)%x(1)=rand4*(-box_size/2.+box_size*real(2*j-1)/(2.*line_size))
          f(line_position)%x(2)=rand1+rand5
          f(line_position)%x(3)=rand2+rand6
          if(j==1) then
            f(line_position)%behind=i*line_size
            f(line_position)%infront=line_position+1
          else if (j==line_size) then
            f(line_position)%behind=line_position-1
            f(line_position)%infront=(i-1)*line_size+1
          else
            f(line_position)%behind=line_position-1
            f(line_position)%infront=line_position+1
          end if
          f(line_position)%u1=0. ; f(line_position)%u2=0. ; f(line_position)%u3=0.
        end do
      end if
    end do
  end subroutine
!*************************************************************************
  subroutine setup_criss_cross_wave
    implicit none
    real :: wave_number, prefactor
    real :: estimate_speed
    real :: amp, random_shift
    real :: xpos, ypos
    real :: dir, axis
    integer :: pcount_required
    integer :: line_size, line_position
    integer :: i, j, k
    !test run.in parameters, if wrong program will exit
    if (line_count==0) then
      call fatal_error('init.mod:criss_cross_wave', &
      'you have not set a value for line_count in run.in')
    end if
    if (periodic_bc.or.periodic_bc_notx.or.periodic_bc_notxy) then
      !work out the number of particles required for single line
      !given the box size specified in run.i
      pcount_required=line_count*nint(box_size/(0.45*delta)) !100% as waves are added
      write(*,*) 'changing size of pcount to fit with box_length and delta'
      write(*,'(a,i7.1)') ' pcount is now:', pcount_required
      deallocate(f) ; pcount=pcount_required ; allocate(f(pcount))
    else
      call fatal_error('init.mod:setup_criss_cross_wavee', &
      'periodic boundary conditions required')
    end if
    write(*,'(a,i3.1,a)') ' drawing', line_count, ' lines in random directions'
    write(*,'(i4.1,a,a,a)') wave_count, ' ',trim(wave_type),' wave pertubations'
    write(*,'(a,i3.1)') ' starting wavenumber ', wave_start
    write(*,'(a,f9.5)') ' starting amplitude', wave_amp*delta
    write(*,'(a,i3.1)') ' wavenumber separation', wave_skip 
    line_size=int(pcount/line_count)
    !START THE LOOP
    do i=1, line_count
      if (line_count==1) then
        xpos=0. ; ypos=0.
      else
        call random_number(xpos) ; call random_number(ypos)
        call random_number(dir) ; call random_number(axis)
        xpos=(xpos-.5)*box_size
        ypos=(ypos-.5)*box_size
        if (dir<0.5) then
          dir=1. 
        else
          dir=-1.
        end if
      end if 
      do j=1, line_size
        line_position=j+(i-1)*line_size
        if (axis<0.333333333) then
          f(line_position)%x(1)=xpos
          f(line_position)%x(2)=ypos
          f(line_position)%x(3)=dir*box_size/2.+box_size*real(2*j-1)/(2.*line_size)
        else if (axis<0.66666666666) then
          f(line_position)%x(2)=xpos
          f(line_position)%x(3)=ypos
          f(line_position)%x(1)=dir*box_size/2.+box_size*real(2*j-1)/(2.*line_size)
        else
          f(line_position)%x(3)=xpos
          f(line_position)%x(1)=ypos
          f(line_position)%x(2)=dir*box_size/2.+box_size*real(2*j-1)/(2.*line_size)
        end if
        if(j==1) then
          f(line_position)%behind=i*line_size
          f(line_position)%infront=line_position+1
        else if (j==line_size) then
          f(line_position)%behind=line_position-1
          f(line_position)%infront=(i-1)*line_size+1
        else
          f(line_position)%behind=line_position-1
          f(line_position)%infront=line_position+1
        end if
        f(line_position)%u1=0. ; f(line_position)%u2=0. ; f(line_position)%u3=0.  
      end do
      !we have now drawn the basic line, now we add the wave pertubations
      prefactor=wave_amp/(wave_start**(wave_slope/2.)) !our starting wavenumber is wave_start
      if (i==1) then !only write this once
        write(*,*)'wave information recorded in ./data/wave_info.log'
        open(unit=34,file='./data/wave_info.log')
        write(34,*) '%------------------WAVE INFORMATION-------------------'
      end if
      do k=1, wave_count !wave_count set in run.in
        wave_number=wave_start+(k-1)*wave_skip
        call random_number(random_shift) !help things along with a 
        random_shift=random_shift*2*pi   !random shift \in (0,2\pi)
        if (wave_count==1) random_shift=0. !0 this for a single wave
        amp=prefactor*(wave_number**(wave_slope/2.))
        if (i==1) then
          write(34,'(f9.5,f9.5,f9.5)') wave_number, amp, random_shift
        end if
        do j=1, line_size !reloop over particles
          line_position=j+(i-1)*line_size
          if (axis<0.333333333333) then
            f(line_position)%x(1)=f(line_position)%x(1)+&
            amp*delta*cos(random_shift+wave_number*2.*pi*real(2.*j-1)/(2.*line_size))
            f(line_position)%x(2)=f(line_position)%x(2)+&
            amp*delta*sin(random_shift+wave_number*2.*pi*real(2.*j-1)/(2.*line_size))
          else if (axis<0.66666666666) then
            f(line_position)%x(2)=f(line_position)%x(2)+&
            amp*delta*cos(random_shift+wave_number*2.*pi*real(2.*j-1)/(2.*line_size))
            f(line_position)%x(3)=f(line_position)%x(3)+&
            amp*delta*sin(random_shift+wave_number*2.*pi*real(2.*j-1)/(2.*line_size))
          else
            f(line_position)%x(3)=f(line_position)%x(3)+&
            amp*delta*cos(random_shift+wave_number*2.*pi*real(2.*j-1)/(2.*line_size))
            f(line_position)%x(1)=f(line_position)%x(1)+&
            amp*delta*sin(random_shift+wave_number*2.*pi*real(2.*j-1)/(2.*line_size))
          end if
        end do
      end do !close k loop
      if (i==1) then
        close(34)
      end if 
    end do !closes the i loop
  end subroutine
  !****************************************************************************
  !>set up a single line from -z to z, number of particles is automatically adjusted
  !>to box size and \f$\delta\f$
  subroutine setup_kursa_recon
    implicit none
    integer :: pcount_required
    integer :: mid_point,mid_point2
    real :: yshift
    real :: angler=0.02
    integer :: i
    if (periodic_bc) then
      !work out the number of particles required for single line
      !given the box size specified in run.i
      pcount_required=2*nint(box_size/(0.75*delta)) !75%
      print*, 'changing size of pcount to fit with box_length and delta'
      print*, 'pcount is now', pcount_required
      deallocate(f) ; pcount=pcount_required ; allocate(f(pcount))
    else
      call fatal_error('init.mod:setup_crow', &
      'periodic boundary conditions required')
    end if
    mid_point=(pcount/2+1)/2
    mid_point2=mid_point+pcount/2
    yshift=angler*0.06*box_size
    write(*,*) 'initf: crow, separation of lines is:', 2.*delta
    !loop over particles setting spatial and 'loop' position
    do i=1, pcount/2
      f(i)%x(1)=0.
      f(i)%x(3)=-box_size/2.+box_size*real(2*i-1)/(pcount)
      if (i<=mid_point) then
        f(i)%x(2)=-yshift+angler*f(i)%x(3)
      else
        f(i)%x(2)=f(mid_point)%x(2)-angler*f(i)%x(3)
      end if  
      if (i==1) then
        f(i)%behind=pcount/2 ; f(i)%infront=i+1
      else if (i==pcount/2) then 
        f(i)%behind=i-1 ; f(i)%infront=1
      else
        f(i)%behind=i-1 ; f(i)%infront=i+1
      end if
      !zero the stored velocities
      f(i)%u1=0. ; f(i)%u2=0.; f(i)%u3=0.
    end do
    !second line
    do i=pcount/2+1, pcount
      f(i)%x(1)=0.
      f(i)%x(3)=box_size/2.-box_size*real(2*(i-pcount/2)-1)/(pcount)
      if (i<=mid_point2) then
        f(i)%x(2)=yshift+angler*f(i)%x(3)
      else
        f(i)%x(2)=f(mid_point2)%x(2)-angler*f(i)%x(3)
      end if  
      if (i==(pcount/2+1)) then
        f(i)%behind=pcount ; f(i)%infront=i+1
      else if (i==pcount) then 
        f(i)%behind=i-1 ; f(i)%infront=1+pcount/2
      else
        f(i)%behind=i-1 ; f(i)%infront=i+1
      end if
      !zero the stored velocities
      f(i)%u1=0. ; f(i)%u2=0. ; f(i)%u3=0.
    end do    
!    print*, 'argggghhhhhhh'
!    print*, f(mid_point)%x(2), f(1)%x(2), angler
!    print*, (f(mid_point)%x(2)-f(1)%x(2))/angler
  end subroutine
  !****************************************************************************
end module
