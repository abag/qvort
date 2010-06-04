module initial
  use cdata
  contains
  !*************************************************************************
  subroutine init_setup()
    implicit none
    pcount=init_pcount
    allocate(f(pcount)) !main vector allocated
    !periodic bounday conditions?
    if (box_size>0.) then
      periodic_bc=.true.
      print*, 'running with periodic boundaries, box size:', box_size
    else 
      print*, 'running with open boundaries'
    end if
    !sort out if there is a special dump time
    int_special_dump=int(special_dump/dt) !convert to integer
    !choose the correct setup routine based on the value of initf in run.in
    select case(initf)
      case('single_loop')
        call setup_single_loop
      case('random_loops')
        call setup_random_loops
      case('leap-frog')
        call setup_leap_frog
      case('line_motion')
        call setup_line_motion
      case('tangle')
        call setup_tangle
      case default
        print*, 'invalid choice for initf parameter'
        print*, 'available options: single_loop, leap-frog'
        print*, 'line_motion, random_loops, tangle'
        call fatal_error !cdata.mod
     end select
  end subroutine
  !*************************************************************************
  subroutine setup_single_loop
    !a loop in the x-y plane
    implicit none
    real :: radius
    integer :: i 
    radius=(0.75*pcount*delta)/(2*pi) !75% of potential size 
    print*, 'initf: single loop, radius of loop:', radius 
    !loop over particles setting spatial and 'loop' position
    do i=1, pcount
      f(i)%x(1)=radius*sin(pi*real(2*i-1)/pcount)
      f(i)%x(2)=radius*cos(pi*real(2*i-1)/pcount) ; f(i)%x(3)=0.
      if (i==1) then
        f(i)%behind=pcount ; f(i)%infront=i+1
      else if (i==pcount) then 
        f(i)%behind=i-1 ; f(i)%infront=1
      else
        f(i)%behind=i-1 ; f(i)%infront=i+1
      end if
      !zero the stored velocities
      f(i)%u1=0. ; f(i)%u2=0.
    end do   
  end subroutine
  !*************************************************************************
  subroutine setup_leap_frog
    !two loops in the x-y plane
    implicit none
    real :: radius
    integer :: i 
    if (mod(pcount,2)/=0) then
      print*, 'pcount is not a multiple of 2-aborting'
      stop
    end if
    radius=(0.75*pcount*delta)/(4*pi) !75% of potential size
    print*, 'initf: leap-frog, radius of loops:', radius 
    !loop over particles setting spatial and 'loop' position
    do i=1, pcount/2
      f(i)%x(1)=radius*sin(pi*real(2*i-1)/(pcount/2))
      f(i)%x(2)=radius*cos(pi*real(2*i-1)/(pcount/2)) ; f(i)%x(3)=0.
      if (i==1) then
        f(i)%behind=pcount/2 ; f(i)%infront=i+1
      else if (i==pcount/2) then 
        f(i)%behind=i-1 ; f(i)%infront=1
      else
        f(i)%behind=i-1 ; f(i)%infront=i+1
      end if
      !zero the stored velocities
      f(i)%u1=0. ; f(i)%u2=0.
    end do
    !second loop
    do i=pcount/2+1, pcount
      f(i)%x(1)=radius*sin(pi*real(2.*(i-pcount/2)-1)/(pcount/2))
      f(i)%x(2)=radius*cos(pi*real(2.*(i-pcount/2)-1)/(pcount/2)) 
      f(i)%x(3)=4.*delta
      if (i==(pcount/2+1)) then
        f(i)%behind=pcount ; f(i)%infront=i+1
      else if (i==pcount) then 
        f(i)%behind=i-1 ; f(i)%infront=1+pcount/2
      else
        f(i)%behind=i-1 ; f(i)%infront=i+1
      end if
      !zero the stored velocities
      f(i)%u1=0. ; f(i)%u2=0.
    end do    
  end subroutine
  !*************************************************************************
  subroutine setup_line_motion
    !a line from the lop of the box to the bottom,
    !given none zero curvature so it moves
    implicit none
    real :: helper_x
    integer :: pcount_required
    integer :: i
    if (periodic_bc) then
      !work out the number of particles required for single line
      !given the box size specified in run.i
      pcount_required=nint(box_size/(0.75*delta)) !75%
      print*, 'changing size of pcount to fit with box_length and delta'
      print*, 'pcount is now', pcount_required
      deallocate(f) ; pcount=pcount_required ; allocate(f(pcount))
    else
      print*, 'periodic boundary conditions required for this initial condition'
      stop
    end if
    do i=1, pcount
      f(i)%x(1)=(box_size/10.)*sin(pi*real(2*i-1)/(2.*pcount))
      f(i)%x(2)=0.
      f(i)%x(3)=-box_size/2.+box_size*real(2*i-1)/(2.*pcount)
      if (i==1) then
        f(i)%behind=pcount ; f(i)%infront=i+1
      else if (i==pcount) then 
        f(i)%behind=i-1 ; f(i)%infront=1
      else
        f(i)%behind=i-1 ; f(i)%infront=i+1
      end if
      !zero the stored velocities
      f(i)%u1=0. ; f(i)%u2=0.
    end do
  end subroutine
  !*************************************************************************
  subroutine setup_random_loops
    implicit none
    real :: loop_radius
    real :: anglex,angley,anglez
    real,dimension(3)::translate
    real,dimension(3)::dummy_xp_1, dummy_xp_2, dummy_xp_3, dummy_xp_4
    integer :: loop_size
    integer:: loop_position
    integer :: i,j
      !test run.in parameters, if wrong program will exit
      if (line_count==0) then
        print*, 'you have not set a value for line_count in run.in'
      end if
      if (mod(pcount,line_count)/=0) then
        print*, 'pcount/line_count is not an integer'
        stop
      end if
      loop_size=int(pcount/line_count)
      loop_radius=loop_size*(0.75*delta)/(2*pi) !75% of potential size
      do i=1, line_count
        call random_number(anglex)
        call random_number(angley)
        call random_number(anglez)
        call random_number(translate)
        anglex=anglex*2*pi
        angley=angley*2*pi
        anglez=anglez*2*pi
        translate=(2*box_size*translate-box_size)
          
        do j=1, loop_size

          loop_position=j+(i-1)*loop_size

          dummy_xp_1(1)=loop_radius*sin(pi*real(2*j-1)/loop_size)
          dummy_xp_1(2)=loop_radius*cos(pi*real(2*j-1)/loop_size)
          dummy_xp_1(3)=0.

          dummy_xp_2(1)=dummy_xp_1(1)
          dummy_xp_2(2)=dummy_xp_1(2)*cos(anglex)+dummy_xp_1(3)*sin(anglex)
          dummy_xp_2(3)=dummy_xp_1(2)*(-sin(anglex))+dummy_xp_1(3)*cos(anglex)

          dummy_xp_3(1)=dummy_xp_2(1)*cos(angley)+dummy_xp_2(3)*(-sin(angley))
          dummy_xp_3(2)=dummy_xp_2(2)
          dummy_xp_3(3)=dummy_xp_2(1)*sin(angley)+dummy_xp_2(3)*cos(angley)

          dummy_xp_4(1)=dummy_xp_3(1)*cos(anglez)+dummy_xp_3(2)*sin(anglez)
          dummy_xp_4(2)=dummy_xp_3(1)*(-sin(anglez))+dummy_xp_3(2)*cos(anglez)
          dummy_xp_4(3)=dummy_xp_3(3)
     
          f(loop_position)%x(1)=dummy_xp_4(1)+translate(1)
          f(loop_position)%x(2)=dummy_xp_4(2)+translate(2)
          f(loop_position)%x(3)=dummy_xp_4(3)+translate(3)

          if(j==1) then
            f(loop_position)%behind=i*loop_size
            f(loop_position)%infront=loop_position+1
          else if (j==loop_size) then
            f(loop_position)%behind=loop_position-1
            f(loop_position)%infront=(i-1)*loop_size+1
          else
            f(loop_position)%behind=loop_position-1
            f(loop_position)%infront=loop_position+1
          end if
          f(loop_position)%u1=0. ; f(loop_position)%u2=0.
        end do
      end do
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
        print*, 'you have not set a value for line_count in run.in'
      end if
      if (mod(pcount,line_count)/=0) then
        print*, 'pcount/line_count is not an integer'
        stop
      end if
      if (periodic_bc) then
        !work out the number of particles required for our lines
        !given the box size specified in run.i
        pcount_required=line_count*nint(box_size/(0.75*delta)) !75%
        print*, 'changing size of pcount to fit with box_length and delta'
        print*, 'pcount is now', pcount_required
        deallocate(f) ; pcount=pcount_required ; allocate(f(pcount))
      else
        print*, 'periodic boundary conditions required for this initial condition'
        call fatal_error !cdata.mod
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
          f(line_position)%u1=0. ; f(line_position)%u2=0.
        end do
      end do
  end subroutine
!****************************************************************
end module
