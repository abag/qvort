module mag
  !MAGNETIC FIELD 
  use cdata
  use general
  use smoothing
  use tree
  contains
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
  subroutine B_smooth()
    !use greens function for diffusion operator to smooth the field
    !use a 5 point footprint
    implicit none
    real, allocatable :: B(:) ! a dummy field to populate with smoothed B
    real, dimension(5) :: dist, greenf !greens function coefficients
    integer :: behind, bbehind, infront, iinfront !helpers
    integer :: i ! for looping
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
        greenf(:)=exp(-(dist(:)**2)/(4*B_nu*dt))
        if (B_3D_nu) then
          greenf=greenf/(3*sum(greenf)-2*greenf(3))
        else
          greenf=greenf/sum(greenf)
        end if
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
