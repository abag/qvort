!> all the routines to employ forcing in the code note this is different to the normal velocity for a vortex system
!> normal fluid enters as \f$\mathbf{u}=\alpha \mathbf{s}' \times (\mathbf{u}_\mathrm{n}-\mathbf{u}_\mathrm{s})+
!>\alpha' \mathbf{s}' \times [\mathbf{s}' \times (\mathbf{u}_\mathrm{n}-\mathbf{u}_\mathrm{s})]\f$ whereas forcing is directly
!>added to vortex velocity i.e.
!>\f$\mathbf{u}=\mathbf{u}_\textrm{vortex}+\mathbf{u}_\mathrm{force}\f$.
!>See \ref FORCE for more details.
module forcing
  use cdata
  use general
  !> @param force_direction a helper vector to force a particle in 3 spatial dimensions
  real, private :: force_direction(3)=0. 
  !> @param force_phase a vector to force Kelvin wave cascade
  real, private :: force_phase(6)=0. 
  real, private :: wave_tail_freq(3)=0.
  !> @param LS_k Large scale forcing wavenumber
  !> @param LS_A Large scale forcing - vector perp to k
  !> @param LS_B Large scale forcing - vector perp to k
  real, private :: LS_k(3), LS_A(3), LS_B(3)
  !> @param LS_k2 Large scale forcing squared wavenumber
  real :: LS_k2
  contains
  !>check all the necessary conditions to use forcing are set in run.in
  subroutine setup_forcing()
    implicit none
    integer :: j !for looping
    select case(force)
      case('off')
        write(*,*) 'no forcing employed'
      case('xflow')
        write(*,*) 'adding an imposed flow in the x direction'
        write(*,'(a,f6.3)') ' u(x)=', force_amp
      case('top_boundary')
        !check initial conditions
         select case(initf)
           case('single_line', 'line_motion', 'lattice')
             write(*,'(a,f10.5,a,f10.5)') 'forcing top boundary with amplitude ', force_amp, &
             ' angular frequency ', force_freq
             write(*,'(a,f13.7)') 'to force at k=10 set frequency to ', &
             (quant_circ*((10.*2.*pi/box_size)**2)/(4*pi))*log((2./(((10.*2.*pi/box_size)**2)*corea))-0.57721)
           case default
             call fatal_error('forcing.mod:setup_forcing', &
             'incorrect initf for forcing to be applied') !cdata.mod
         end select
      case('box_shake')
        write(*,'(a,f5.3,a,i4.3,a)') 'box shaking forcing with amplitude ', force_amp, &
        ' forcing direction changed every ', ceiling(force_freq), 'timesteps'
      case('potential_vortex')
        write(*,'(a)') ' simple potential vortex in xy-plane'
      case('looping_potential_vortex')
        write(*,'(a)') ' potential vortex in xy-plane that moves in a loop'
      case('delta_corr')
        write(*,'(a,f5.3)') 'delta correlated (time/space) forcing with amplitude ', force_amp
      case('LS_force')
        write(*,'(a,f5.3)') 'large scale forcing (similar to KS) with amplitude ', force_amp
      case('wave_force')
        write(*,'(a,f5.3)') 'force wavemodes  9,10,11 with amplitude ', force_amp
        write(*,'(a,i6.2,a)') 'phase changed every, ', ceiling(force_freq), ',timesteps'
      case('wave_tail')
        write(*,'(a,f5.3)') 'force tail of vortex  with modes 9,10,11 with amplitude ', force_amp
        write(*,'(a,i6.2,a)') 'phase changed every, ', ceiling(force_freq), ',timesteps'
        do j=9,11
          wave_tail_freq(j-8)=(quant_circ*((j*2.*pi/box_size)**2)/(4*pi))*&
                                   log((2./(((10.*2.*pi/box_size)**2)*corea))-0.57721)
        end do
        write(*,'(a,3f13.6)') 'frequencies are, ', wave_tail_freq(:)
      case default
        call fatal_error('forcing.mod:setup_forcing', &
        'incorrect forcing parameter set in run.in') !cdata.mod
      end select
      if (force_cutoff<1000.) then
        write(*,'(a,f6.3)') ' forcing is turned off when t=',force_cutoff 
      end if
      if (forcing_mesh_fileprint) write(*,*) 'printing forcing array to file'
  end subroutine
  !***********************************************
  !>force the particles - an additional velocity at position
  !>of f(i)%x
  subroutine get_forcing(i,u,u_BS)
    implicit none
    integer, intent(IN) :: i
    real, intent(OUT) :: u(3)
    real, intent(INOUT) :: u_BS(3)
    integer :: j !for looping
    integer :: peri, perj, perk !used to loop in periodic cases
    select case(force)
      case('off')
        u=0.
      case('xflow')
        u=0.
        u(1)=force_amp
      case('top_boundary')
        u=0.
        if ((f(i)%x(3)-box_size/2.)>-1.5*delta) then
          !particle is sufficiently close to top boundary to force
          u=0.
          !sinusoidal forcing in the x direction
          !in run.in we are setting angular frequency
          u(1)=force_amp*sin(force_freq*t)
        end if
      case('box_shake')
        u=force_direction*force_amp
      case('potential_vortex')
        u(1)=-f(i)%x(2)/(f(i)%x(1)**2+f(i)%x(2)**2+(box_size/80)**2)
        u(2)=f(i)%x(1)/(f(i)%x(1)**2+f(i)%x(2)**2+(box_size/80)**2)
        u(3)=0.
        if (periodic_bc) then
          do peri=-1,1 ; do perj=-1,1 
            if (peri==0.and.perj==0) cycle
            u(1)=u(1)-(f(i)%x(2)-perj*box_size)/&
            ((f(i)%x(1)-peri*box_size)**2+(f(i)%x(2)-perj*box_size)**2+(box_size/80)**2)
            u(2)=u(2)+(f(i)%x(1)-peri*box_size)/&
            ((f(i)%x(1)-peri*box_size)**2+(f(i)%x(2)-perj*box_size)**2+(box_size/80)**2)
          end do ; end do
        end if
        u=u*(box_size/40)
      case('looping_potential_vortex')
        u(1)=-(f(i)%x(2)-force_direction(2))/&
            ((f(i)%x(1)-force_direction(1))**2+(f(i)%x(2)-force_direction(2))**2+(box_size/80)**2)
        u(2)=(f(i)%x(1)-force_direction(1))/&
            ((f(i)%x(1)-force_direction(1))**2+(f(i)%x(2)-force_direction(2))**2+(box_size/80)**2)
        u(3)=0.
        if (periodic_bc) then
          do peri=-1,1 ; do perj=-1,1 
            if (peri==0.and.perj==0) cycle
            u(1)=u(1)-(f(i)%x(2)-perj*box_size-force_direction(2))/&
            ((f(i)%x(1)-peri*box_size-force_direction(1))**2+&
            (f(i)%x(2)-perj*box_size-force_direction(2))**2+(box_size/80)**2)
            u(2)=u(2)+(f(i)%x(1)-peri*box_size-force_direction(1))/&
            ((f(i)%x(1)-peri*box_size-force_direction(1))**2+&
            (f(i)%x(2)-perj*box_size-force_direction(2))**2+(box_size/80)**2)
          end do ; end do
        end if
        u=u*(box_size/40)
      case('delta_corr')
        !at each time and position generate a new random vector
        call random_number(force_direction)
        force_direction=force_direction*2.-1.
        !normalise it
        force_direction=force_direction/sqrt(dot_product(force_direction,force_direction))
        !use this as the forcing
        u=force_direction*force_amp
      case('LS_force')
        u=cross_product(LS_A,LS_k)*sin(dot_product(LS_k,f(i)%x))/LS_k2+cross_product(LS_B,LS_k)*cos(dot_product(LS_k,f(i)%x))/LS_k2
        !use this as the forcing-multiply by amplitude
        u=u*force_amp
      case('wave_force')
        u=0. !initialise to 0
        !now sum over three wavenumbers
        do j=9,11
          u(1)=u(1)+force_amp*cos(2*pi*j*f(i)%x(3)/box_size+force_phase(j-8))
          u(2)=u(2)+force_amp*sin(2*pi*j*f(i)%x(3)/box_size+force_phase(j-8))
          u(1)=u(1)+force_amp*cos(-2*pi*j*f(i)%x(3)/box_size+force_phase(j-5))
          u(2)=u(2)+force_amp*sin(-2*pi*j*f(i)%x(3)/box_size+force_phase(j-5))
        end do
      case('wave_tail')
        u=0.
        if (abs(f(i)%x(3))<epsilon(0.)) then
          !print*, i
          !0 the Biot-Savart velocity
          u_BS=0.
          !particle is sufficiently close to top boundary to force
          open(unit=67,file='./data/wave_tail_pos.log',position='append')
            write(67,*) i,itime, f(i)%x(1), f(i)%x(2)
          close(67)
          do j=1,1
            !u(1)=u(1)+force_amp*cos(t*wave_tail_freq(j))
            u(2)=u(2)+force_amp*sin(t*wave_tail_freq(j)+pi/2)
          end do
        end if
    end select
  end subroutine
  !***********************************************
  !>routine to randomise certain forcing routines every timestep, if required
  !>we can also print to file in here
  subroutine randomise_forcing
    implicit none
    integer :: j
    select case(force)
      case('box_shake')
        if (mod(itime,ceiling(force_freq))==0) then
          call random_number(force_direction)
          force_direction=force_direction*2.-1.
          !normalise
          force_direction=force_direction/sqrt(dot_product(force_direction,force_direction))
        end if
      case('LS_force')
        !this has it's own subroutine to reinitialise
        if (mod(itime,ceiling(force_freq))==0) then
          call get_LS_kAB!forcing.mod
          if (mod(itime,shots)==0) then
            !print to file
            open(unit=47,file='./data/LS_forcing.log',position='append')
              write(47,*) LS_k, LS_A, LS_B
            close(47)
          end if 
        end if
      case('looping_potential_vortex')
        !the potential vortex 
        force_direction(1)=(box_size/4.)*cos(2*pi*t*force_freq)
        force_direction(2)=(box_size/4.)*sin(2*pi*t*force_freq)
        force_direction(3)=0.
        if (mod(itime,shots)==0) then
           open(unit=72,file='./data/looping_pot_vortex.log',position='append')
            write(72,*) t, force_direction(1:3)
          close(72)
        end if
      case('wave_force','wave_tail')
        !uniform distributed random phase from 0->2\pi
        if (mod(itime,ceiling(force_freq))==0) then
          do j=1, 6
            force_phase(j)=runif(0.,2.*pi)
          end do 
        end if
    end select 
    if (forcing_mesh_fileprint) then
      if (mod(itime,mesh_shots)==0) then
        call print_forcing_mesh
      end if
    end if
  end subroutine
    !**********************************************************
    !>if selected print forcing field (discretized on a cubic mesh) to a binary file
    subroutine print_forcing_mesh
      implicit none
      real :: forcing_mesh(64,64,3) !only 2d at present
      real :: x(64), u(3)
      integer :: i, j !for looping
      character (len=40) :: print_file
      do i=1, 64
        x(i)=(box_size/64)*real(2*i-1)/2.-(box_size/2.)
      end do 
      do j=1, 64  ; do i=1, 64
        !get the forcing field 
        call get_forcing_gen((/x(i),x(j),0./),u) !forcing must be in general forcing
        forcing_mesh(j,i,1)=u(1)
        forcing_mesh(j,i,2)=u(2)
        forcing_mesh(j,i,3)=u(3)
      end do ; end do 
      write(unit=print_file,fmt="(a,i3.3,a)")"./data/forcing_mesh",itime/mesh_shots,".dat"
      open(unit=92,file=print_file,form='unformatted',status='replace',access='stream')
        write(92) x(1:64)
        write(92) forcing_mesh(:,:,1)
        write(92) forcing_mesh(:,:,2)
        write(92) forcing_mesh(:,:,2)
      close(92)
    end subroutine
  !***********************************************
  !>Large scale forcing routine, very similar to KS with one
  !>random large scale mode, this sets up k, A and B for use in
  !>main forcing routine
  subroutine get_LS_kAB
    real,dimension(3) ::  newa, newb, unitk !helper 
    integer :: i
    !create random vectors
    call random_number(newa) ; call random_number(newb)
    newa=2.*newa-1. ; newa2=2.*newa2-1.
    newa=newa/(sqrt(sum(newa**2))) ; newb=newb/(sqrt(sum(newb**2)))
    !generate random wavevector
    LS_k=.0
    do while (sum(LS_k)<epsilon(0.))
      do i=1,3
        call random_number(LS_k(i))
        if (LS_k(i)<0.33333) then
          LS_k(i)=1.
        else if (LS_k(i)<0.666) then
          LS_k(i)=-1
        else
          LS_k(i)=0.
        end if
      end do
    end do
    LS_k=LS_k*2*pi/box_size !scale by 2\pi for periodicity - box size must enter here
    LS_k2=sqrt(sum(LS_k**2))
    unitk=LS_k/LS_k2 !unit vector
    !create A and B by taking cross product of newa/b with unit k
    LS_A=cross_product(newa,unitk) ; LS_B=cross_product(newb,unitk)
  end subroutine
  !***********************************************
  !>a general forcing routine which accepts a position
  subroutine get_forcing_gen(x,u)
    implicit none
    real, intent(IN) :: x(3)
    real, intent(OUT) :: u(3)
    integer :: peri, perj, perk !used to loop in periodic cases
    select case(force)
      case('off')
        u=0.
      case('top_boundary')
        if ((x(3)-box_size/2.)>-1.5*delta) then
          !particle is sufficiently close to top boundary to force
          u=0.
          !sinusoidal forcing in the x direction
          u(1)=delta*force_amp*sin(force_freq*t/(2*pi))
        end if
      case('delta_corr')
        !at each time and position generate a new random vector
        call random_number(force_direction)
        force_direction=force_direction*2.-1.
        !normalise it
        force_direction=force_direction/sqrt(dot_product(force_direction,force_direction))
        !use this as the forcing
        u=force_direction*force_amp
      case('potential_vortex')
        u(1)=-x(2)/(x(1)**2+x(2)**2+(box_size/80)**2)
        u(2)=x(1)/(x(1)**2+x(2)**2+(box_size/80)**2)
        u(3)=0.
        if (periodic_bc) then
          do peri=-1,1 ; do perj=-1,1 
            if (peri==0.and.perj==0) cycle
            u(1)=u(1)-(x(2)-perj*box_size)/&
            ((x(1)-peri*box_size)**2+(x(2)-perj*box_size)**2+(box_size/80)**2)
            u(2)=u(2)+(x(1)-peri*box_size)/&
            ((x(1)-peri*box_size)**2+(x(2)-perj*box_size)**2+(box_size/80)**2)
          end do ; end do
        end if
        u=u*(box_size/40)
      case('looping_potential_vortex')
        u(1)=-(x(2)-force_direction(2))/&
            ((x(1)-force_direction(1))**2+(x(2)-force_direction(2))**2+(box_size/80)**2)
        u(2)=(x(1)-force_direction(1))/&
            ((x(1)-force_direction(1))**2+(x(2)-force_direction(2))**2+(box_size/80)**2)
        u(3)=0.
        if (periodic_bc) then
          do peri=-1,1 ; do perj=-1,1 
            if (peri==0.and.perj==0) cycle
            u(1)=u(1)-(x(2)-perj*box_size-force_direction(2))/&
            ((x(1)-peri*box_size-force_direction(1))**2+&
            (x(2)-perj*box_size-force_direction(2))**2+(box_size/80)**2)
            u(2)=u(2)+(x(1)-peri*box_size-force_direction(1))/&
            ((x(1)-peri*box_size-force_direction(1))**2+&
            (x(2)-perj*box_size-force_direction(2))**2+(box_size/80)**2)
          end do ; end do
        end if
        u=u*(box_size/40)
      case default
        call fatal_error('forcing.mod','this forcing function is not compatible with a general call')
    end select
  end subroutine
end module
!> \page FORCE Forcing
!!Forcing of filaments is set in run.in throught the parameter
!!force this accepts the following arguements\n
!!
!!\p zero - no forcing, this is the default \n
!!\p xflow - imposed flow in the x direction, velocity set by force_amp \n
!!\p box_shake - correlated forcing moving box randomly with frequency set by
!!\p force_feq in run.in \n
!!\p delta_corr - delta correlated forcing in time and space with amplitude set
!!by \p force_amp parameter \n
!!\p top_boundary - sinusoidal forcing in the x direction at the top of the box
!!with a frequency and amplitude set in run.in (force_freq, force_amp)\n
!!\p LS_force - please document me\n
!!\p wave_force - forces at three specified wavenumbers
!!\f[ 
!!u_{x,{\rm force}}(z)=\sum_{k=9}^{11} \Re [ A \exp{i(kz+\phi_i)} ] =\sum_{k=9}^{11} A\cos(kz+\phi_i), \\
!!u_{y,{\rm force}}(z)=\sum_{k=9}^{11} \Im [ A \exp{i(kz+\phi_i)} ] =\sum_{k=9}^{11} A\sin(kz+\phi_i),
!!\f] 
!!where \f$ \phi\f$ is a randomly chosen phase, reset every force_freq timesteps.

