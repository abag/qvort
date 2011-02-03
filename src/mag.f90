module mag
  !MAGNETIC FIELD 
  use cdata
  contains
  !**********************************************************************
  subroutine B_ts
    implicit none
    open(unit=71,file='data/B_ts.log',position='append')
      if (itime==shots) then
        write(71,*) '----t-------Bmax-------Bmin----'
      end if
      write(71,'(f13.7,f13.7,f13.7)') t, maxval(f(:)%B), minval(f(:)%B)
    close(71)
  end subroutine
end module
