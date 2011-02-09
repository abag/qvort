module mag
  !MAGNETIC FIELD 
  use cdata
  use general
  contains
  !**********************************************************************
  subroutine mag_tension(i,u)
    implicit none
    integer,intent(IN) :: i !particle
    real, intent(OUT) :: u(3) !tension velocity
    call get_deriv_2(i,u)
    u=u*(f(i)%B**2)
  end subroutine
  !**********************************************************************
  subroutine B_ts
    implicit none
    open(unit=71,file='data/B_ts.log',position='append')
      if (itime==shots) then
        write(71,*) '%--------t--------Bmax---------Bmin------------Brms-------'
      end if
      write(71,'(f13.7,f13.7,f13.7,f13.7)') t, maxval(f(:)%B), minval(f(:)%B), Brms
    close(71)
  end subroutine
  !**********************************************************************
  subroutine Bstretch(filenumber)
    implicit none
    integer, intent(IN) :: filenumber
    character (len=40) :: print_file
    integer :: i
    write(unit=print_file,fmt="(a,i4.4,a)")"./data/Bstretch",filenumber,".log"
    open(unit=98,file=print_file,status='replace')
    do i=1,pcount
      if (f(i)%infront==0) cycle
      write(98,*) f(i)%B
    end do
    close(98)
  end subroutine
end module
