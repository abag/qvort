!>all routines related to the filaments acting as magnetic flux tubes are in this
!>module, see \ref MAG for details of the algorithms used.
module mag
  use cdata
  use general
  use smoothing
  use tree
  use normal_fluid
  !>@param diffusion_on a flag set on if diffusion had been set to a none
  !>zero value in run.in, this helps subsequent smoothing routine
  logical, private :: diffusion_on=.false.
  real, private :: B_E=0. !magnetic energy
  contains
  !**********************************************************************
  !>sets up the magnetic field and check there are no conflicting options
  !>arguments set in run.in
  subroutine setup_mag()
    implicit none
    write(*,'(a)') ' ---------------------ACTING AS A MAGNETIC FIELD----------------------'
    write(*,'(a, f6.3, a)') ' initial field strength ', B_init, ' G'
    if (B_nu>epsilon(0.)) then
      diffusion_on=.true.
      write(*,'(a, f10.6, a)') ' magnetic diffusivity ', B_nu, ' cm^3/s'
      !check our magnetic diffusivity cfl condition
      write(*,'(a,f10.4)') ' magnetic CFL number:', delta**2/(B_nu*dt)
      if (B_3D_nu) then
        write(*,*)'diffusion is three dimensional'
      else
        write(*,*)'diffusion only acts along flux ropes'
      end if
    else
      write(*,*) 'magnetic diffusivity is not acting on flux ropes'
    end if
    f(:)%B=B_init
    !"f(1:20)%B=4*B_init
    if (full_B_print) write(*,*) 'printing full B field to file'
    write(*,'(a,e10.4)') ' magnetic tension coefficient ', B_tension
  end subroutine
  !**********************************************************************
  !>get the length between the particle and the particle infront-assign to l1
  subroutine set_B_length1
    implicit none
    integer :: i
    do i=1, pcount
      if (f(i)%infront==0) cycle !empty particle
      !set l1 the distance between the particle and the one infront
      f(i)%l1=dist_gen(f(i)%x,f(i)%ghosti) !general.f90
      f(i)%v1=1. !set this to 1
    end do
  end subroutine
  !**********************************************************************
  !>get the length between the particle and the particle infront-assign to l2
  subroutine set_B_strength
    use sph_interface
    implicit none
    real :: divu
    integer :: i
    B_E=0.
    do i=1, pcount
      if (f(i)%infront==0) cycle !empty particle
      !set l2 the distance between the particle and the one infront
      f(i)%l2=dist_gen(f(i)%x,f(i)%ghosti) !general.f90
      !account for compressibility of velocity field
      if (nf_compressible) then
        !get the local divergence of the velocity field
        select case(velocity) ; case('SPH')
          call SPH_f_divu(i,divu) !sph_interface.mod
        case default
          call get_normal_divu(f(i)%x,divu) !normal.mod
        end select
        !now timestep the volume element
        f(i)%v2=f(i)%v1*(1.+dt*divu)
        !now set field strength
        f(i)%B=f(i)%B*(f(i)%l2/f(i)%l1)*(f(i)%v1/f(i)%v2)
        if (f(i)%B>10) then
          print*, i, f(i)%B, (f(i)%l2/f(i)%l1), f(i)%v2/f(i)%v1
          print*, i, divu
        end if
      else
        !now set field strength
        f(i)%B=f(i)%B*(f(i)%l2/f(i)%l1)
      end if
      B_E=B_E+f(i)%B*f(i)%l2
    end do
    select case(normal_velocity)
      case('shear')
        if ((mod(itime,shots)==0)) then
          open(unit=23,file='data/shear_y_test.log', position='append')
            write(23,*) t, f(10)%x(3), f(10)%B, &
                           f(20)%x(3), f(20)%B, &
                           f(44)%x(3), f(44)%B
          close(23)
        end if
    end select
  end subroutine
  !**********************************************************************
  !>simulate the effect of the lorentz force by adding in magnetic tension 
  !>to velocity calculations
  !>\f[\mathbf{J} \times \mathbf{B}=(\nabla \times \mathbf{B}) \times \mathbf{B}
  !>=B^2 d\hat{\mathbf{B}}/ds\f]
  subroutine mag_tension(i,u)
    implicit none
    integer,intent(IN) :: i !particle
    real, intent(OUT) :: u(3) !tension velocity
    real, dimension(3) :: B,Bpx,Bmx,Bpy,Bmy,Bpz,Bmz,J
    logical, parameter :: fullJxB=.false.
    if (fullJxB) then
      B=0. ; call tree_smooth(f(i)%x,vtree,(/0.,0.,0./),B)
      Bpx=0. ; call tree_smooth(f(i)%x+(/delta,0.,0./),vtree,(/0.,0.,0./),Bpx)
      Bmx=0. ; call tree_smooth(f(i)%x-(/delta,0.,0./),vtree,(/0.,0.,0./),Bmx)
      Bpy=0. ; call tree_smooth(f(i)%x+(/0.,delta,0./),vtree,(/0.,0.,0./),Bpy)
      Bmy=0. ; call tree_smooth(f(i)%x-(/0.,delta,0./),vtree,(/0.,0.,0./),Bmy)
      Bpz=0. ; call tree_smooth(f(i)%x+(/0.,0.,delta/),vtree,(/0.,0.,0./),Bpz)
      Bmz=0. ; call tree_smooth(f(i)%x-(/0.,0.,delta/),vtree,(/0.,0.,0./),Bmz)
      J(1)=Bpy(3)-Bmy(3)-Bpz(2)+Bmz(2)
      J(2)=Bpz(1)-Bmz(1)-Bpx(3)+Bmx(3)
      J(3)=Bpx(2)-Bmx(2)-Bpy(1)+Bmy(1)
      J=J/(2*delta)
      u=cross_product(J,B)
    else
      call get_deriv_2(i,u)
      u=u*(f(i)%B**2)
    end if
  end subroutine
    !**********************************************************************
  !>simulate the effect of the lorentz force by adding in magnetic tension 
  !>to velocity calculations, this includes the mirror force:
  !>\f[\mathbf{J} \times \mathbf{B}=B d{\mathbf{B}}/ds\f]
  subroutine mag_tension_mirr(i,u)
    implicit none
    integer,intent(IN) :: i !particle
    real, intent(OUT) :: u(3) !tension velocity
    real, dimension(3) :: B_infront, B_behind
    real :: disti, distb
    disti=dist_gen(f(i)%ghosti,f(i)%x)
    distb=dist_gen(f(i)%ghostb,f(i)%x)
    B_infront=f(i)%B*(f(i)%ghosti-f(i)%x)/disti
    B_behind=f(f(i)%behind)%B*(f(i)%x-f(i)%ghostb)/distb
    u=2.*(B_infront-B_behind)/(disti+distb)
  end subroutine
  !**********************************************************************
  !>use greens function for diffusion operator to smooth the field
  !>use a 5 point footprint, can be done in one or three dimensions
  !>in one dimension this is given by
  !>\f[B(\mathbf{s}_i,t+\Delta t)=\sum_{j=i-2}^{i+2}\exp(-(\mathbf{s}_i
  !>-\mathbf{s}_j)^2/4\nu \Delta t)\frac{B_j}{2\sqrt{\pi\nu\Delta t/\ell_j^2}}\f]
  !>in three dimensions the 2 in the denominator is an 8 select this in run.in
  !>via the parameter B_3d_nu, the magnitude of the diffusion is set by the
  !>parameter B_nu in run.in 
  subroutine B_diffusion()
    implicit none
    real, allocatable :: B(:) ! a dummy field to populate with smoothed B
    real, dimension(5) :: prefact, dist, greenf !greens function coefficients
    integer :: behind, bbehind, infront, iinfront !helpers
    integer :: i, j ! for looping
    if (diffusion_on.eqv..false.) return
    B_E=0. !0 the magnetic energy
    allocate(B(pcount))
   
    do i=1, pcount
      if (f(i)%infront==0) then
        B(i)=0. !empty particle
      else
        behind=f(i)%behind ; bbehind=f(behind)%behind
        infront=f(i)%infront ; iinfront=f(infront)%infront
        dist(2)=dist_gen(f(i)%ghostb,f(i)%x)
        dist(4)=dist_gen(f(i)%ghosti,f(i)%x)
        dist(1)=dist_gen(f(behind)%ghostb,f(behind)%x)+dist(2)
        dist(5)=dist_gen(f(infront)%ghosti,f(infront)%x)+dist(4)
        dist(3)=0. !the particle itself 
        greenf(:)=exp(-(dist(:)**2)/(4*B_nu*dt)) !exponential part of G(x,x')
        prefact(1)=dist(1)-dist(2) !dist_gen(f(bbehind)%x,f(bbehind)%ghosti)
        prefact(2)=dist(2) !dist_gen(f(behind)%x,f(behind)%ghosti)
        prefact(3)=dist(4) !dist_gen(f(i)%x,f(i)%ghosti)
        prefact(4)=dist(5)-dist(4) !dist_gen(f(infront)%x,f(infront)%ghosti)
        prefact(5)=dist_gen(f(iinfront)%x,f(iinfront)%ghosti) !not calc. previously
        prefact=0.5*prefact/sqrt(pi*B_nu*dt) !1D correction
        greenf=greenf*prefact
        if (B_3D_nu) then
          greenf=greenf/(3.*sum(greenf)-2*greenf(3)) !account for 3d diffusion
        else
          greenf=greenf/sum(greenf) !normalise the integral should be 1
        end if
        if (any(isnan(greenf))) then
          print*, itime
          print*, i, behind, infront, bbehind, iinfront
          print*, f(i)%x
          print*, 'ghost infront',f(i)%ghosti!, f(i)%ghostb
          print*, 'ghost behind',f(i)%ghostb
          print*, 'f(infront)%ghosti',f(infront)%ghosti 
          print*, 'dist', dist
          print*, 'greenf',greenf
          call fatal_error('mag.mod:B_diffusion',&
          'serious problem, is B_nu 0?')
        end if
        if (greenf(3)>1.) call fatal_error('mag.mod:B_diffusion',&
          'B_nu is too small for current particle resolution')
        B(i)=greenf(1)*f(bbehind)%B+greenf(2)*f(behind)%B+greenf(3)*f(i)%B+&
             greenf(4)*f(infront)%B+greenf(5)*f(iinfront)%B
      end if
      B_E=B_E+B(i)*dist(4)
    end do
    f(:)%B=B(:)
    !print*, greenf, sum(greenf)
    deallocate(B)
  end subroutine
  !**********************************************************************
  !>print time series information for the magnetic field
  subroutine B_ts
    implicit none
    open(unit=71,file='data/B_ts.log',position='append')
      if (itime==shots) then
        write(71,*) '%--------t--------Bmax---------Bmin------------Brms--------B_energy-----'
      end if
      write(71,'(f13.7,f13.7,f13.7,f13.7, f13.7)') t, maxval(f(:)%B,mask=f(:)%infront>0), minval(f(:)%B,mask=f(:)%infront>0), &
                                            Brms, B_E
      write(*,'(f13.7,f13.7,f13.7,f13.7)',advance='yes')  maxval(f(:)%B,mask=f(:)%infront>0), &
                                                           minval(f(:)%B,mask=f(:)%infront>0), Brms, B_E
    close(71)
  end subroutine
  !**********************************************************************
  !>print full magnetic field information at a particular timestep
  subroutine print_full_B(filenumber)
    implicit none
    integer, intent(IN) :: filenumber
    character (len=40) :: print_file
    integer :: i
    write(unit=print_file,fmt="(a,i4.4,a)")"./data/full_B",filenumber,".dat"
    open(unit=98,file=print_file,status='replace',form='unformatted',access='stream')
      write(98) pcount
      write(98) t
      write(98) f(:)%B
    close(98)
  end subroutine
end module
!>\page MAG magnetic fields
!>Please fill me in!
