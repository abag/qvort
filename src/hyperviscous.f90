!>Routines to add hyperviscosity - enters as and additional velocity
!>\f$\mathbf{u}_{\rm hyp}=\alpha_{hyp} \mathbf{s}' \times -\mathbf{u}_\mathrm{s}\f$.
!>Here \f$  \alpha_{hyp}\f$ is dependent on curvature i.e. \f$ \alpha_{hyp} \propto \kappa^n \f$
module hyperviscous
  use cdata
  use general 
    !count number of times we apply hyperviscosity
    real, private :: hyp_const
    real,private  :: max_curv
    contains 
    !************************************************************
    !>setup everything needed to use hyperviscosity
    !> construct form of \f$ \alpha_{hyp} \f$
    subroutine setup_hyperviscosity
      !in here we print to file the hyperviscous profile
      implicit none
      real :: plot_x(100), plot_y(100)
      integer :: i
      write(*,*) 'hyperviscosity swithed on'
      max_curv=2./delta
      write(*,'(a,f9.2,a)') ' curvature above ', max_curv*hyp_curv, ' are damped'
      write(*,'(a,i4.2)') ' degree of polynomial ', hyp_power
      !a constant to make sure whatever the polynomal alpha(curv_max)=hyp_nu
      hyp_const=hyp_nu/((max_curv-hyp_curv*max_curv)**hyp_power)
      write(*,'(a,e13.6)') ' hyperviscosity: ', hyp_const
      open(unit=37,file='./data/hyperviscous.log')
      do i=1,100
        !create a vector (length 100) from 0 to max curvature
        plot_x(i)=max_curv*(2*i-1)/200.
        if (plot_x(i)<hyp_curv*max_curv) then
          plot_y(i)=0.
        else
          plot_y(i)=hyp_const*(plot_x(i)-hyp_curv*max_curv)**hyp_power
        end if
        !print function to file
        write(37,*) plot_x(i),plot_y(i)
      end do 
      close(37)
    end subroutine
    !************************************************************
    !> Find \f$ \alpha_{hyp} \f$ to return to timestep to implement
    !>hyperviscosity.
    subroutine get_hyp_alpha(curv,hyp_alpha)
      implicit none
      real, intent(IN) :: curv !curvature
      real, intent(OUT) :: hyp_alpha !alpha associated with curvature
      if (curv<hyp_curv*max_curv) then
        hyp_alpha=0.
      else
        !simple form for hyperviscous function defined below
        hyp_alpha=hyp_const*(curv-hyp_curv*max_curv)**hyp_power
      end if
    end subroutine
end module
