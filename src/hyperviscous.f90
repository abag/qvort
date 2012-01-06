!>Routines to add hyperviscosity - complete me!
module hyperviscous
  use cdata
  use general 
    !count number of times we apply hyperviscosity
    integer, private :: hyperviscous_count=0
    real, private :: hyp_const
    real,private  :: max_curv
    contains 
    !************************************************************
    !>setup everything needed to use hyperviscosity
    subroutine setup_hyperviscosity
      !in here we print to file the hyperviscous profile
      implicit none
      real :: plot_x(100), plot_y(100)
      integer :: i
      write(*,*) 'hyperviscosity swithed on'
      max_curv=2./delta
      write(*,'(a,f9.2,a)') ' curvature above ', max_curv*hyp_curv, ' are damped'
      write(*,'(a,i4.2)') ' degree of polynomial ', hyp_power
      hyp_const=1./((max_curv-hyp_curv*max_curv)**hyp_power)
      open(unit=37,file='./data/hyperviscous.log')
      do i=1,100
        plot_x(i)=max_curv*(2*i-1)/200.
        if (plot_x(i)<hyp_curv*max_curv) then
          plot_y(i)=0.
        else
          plot_y(i)=hyp_const*(plot_x(i)-hyp_curv*max_curv)**hyp_power
        end if
        write(37,*) plot_x(i),plot_y(i)
      end do 
      close(37)
    end subroutine
    !************************************************************
    !>
    subroutine get_hyp_alpha(curv,hyp_alpha)
      implicit none
      real, intent(IN) :: curv !curvature
      real, intent(OUT) :: hyp_alpha !alpha associated with curvature
      if (curv<hyp_curv*max_curv) then
        hyp_alpha=0.
      else
        hyp_alpha=hyp_const*(curv-hyp_curv*max_curv)**hyp_power
      end if
    end subroutine
end module
