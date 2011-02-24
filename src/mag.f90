module mag
  !MAGNETIC FIELD 
  use cdata
  use general
  use smoothing
  use tree
  logical, private :: diffusion_on=.false.
  contains
  !**********************************************************************
  subroutine setup_mag()
    implicit none
    write(*,'(a)') ' ---------------------ACTING AS A MAGNETIC FIELD----------------------'
    write(*,'(a, f6.3, a)') ' initial field strength ', B_init, ' G'
    if (B_nu>epsilon(0.)) then
      diffusion_on=.true.
      write(*,'(a, f10.6, a)') ' magnetic diffusivity ', B_nu, ' cm^3/s'
      !check our magnetic diffusivity cfl condition
      write(*,'(a,f10.4)') ' angetic CFL number:', delta**2/(B_nu*dt)
      if (B_3D_nu) then
        write(*,*)'diffusion is three dimensional'
      else
        write(*,*)'diffusion only acts along flux ropes'
      end if
    else
      write(*,*) 'magnetic diffusivity is not acting on flux ropes'
    end if
    f(:)%B=B_init
    if (full_B_print) write(*,*) 'printing full B field to file'
  end subroutine
  !**********************************************************************
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
  subroutine B_diffusion()
    !use greens function for diffusion operator to smooth the field
    !use a 5 point footprint
    implicit none
    real, allocatable :: B(:) ! a dummy field to populate with smoothed B
    real, dimension(5) :: prefact, dist, greenf !greens function coefficients
    integer :: behind, bbehind, infront, iinfront !helpers
    integer :: i ! for looping
    if (diffusion_on.eqv..false.) return
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
        if (B_3D_nu) then
          prefact=0.125*prefact/sqrt(pi*B_nu*dt) !3D correction to above
        else
          prefact=0.5*prefact/sqrt(pi*B_nu*dt) !1D correction
        end if
        greenf=greenf*prefact
        if (any(isnan(greenf))) call fatal_error('mag.mod:B_diffusion',&
        'serious problem, is B_nu 0?')
        if (greenf(3)>1.) call fatal_error('mag.mod:B_diffusion',&
        'B_nu is too small for current particle resolution')
        B(i)=greenf(1)*f(bbehind)%B+greenf(2)*f(behind)%B+greenf(3)*f(i)%B+&
             greenf(4)*f(infront)%B+greenf(5)*f(iinfront)%B
      end if
    end do
    f(:)%B=B(:)
    deallocate(B)
  end subroutine
  !**********************************************************************
  subroutine B_ts
    implicit none
    open(unit=71,file='data/B_ts.log',position='append')
      if (itime==shots) then
        write(71,*) '%--------t--------Bmax---------Bmin------------Brms-------'
      end if
      write(71,'(f13.7,f13.7,f13.7,f13.7)') t, maxval(f(:)%B,mask=f(:)%infront>0), minval(f(:)%B,mask=f(:)%infront>0), Brms
    close(71)
  end subroutine
  !**********************************************************************
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
