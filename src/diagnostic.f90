module diagnostic
  !KEY DIAGNOSTICS - MAINLY TO PRINT TO SCREEN
  !ENERGY REALLY IS INVALID WITH PERIODIC BOUNDARY CONDITIONS
  use cdata
  use general
  contains
  !*************************************************
  subroutine velocity_info()
    implicit none
    real, allocatable :: uinfo(:,:)
    allocate(uinfo(pcount,2))
    uinfo(:,1)=sqrt(f(:)%u(1)**2+f(:)%u(2)**2+f(:)%u(3)**2)
    uinfo(:,2)=sqrt((f(:)%u1(1)-f(:)%u2(1))**2+&
                    (f(:)%u1(2)-f(:)%u2(2))**2+&
                    (f(:)%u1(3)-f(:)%u2(3))**2)
    maxu=maxval(uinfo(:,1)) ; maxdu=maxval(uinfo(:,2))
    deallocate(uinfo)
  end subroutine
  !*************************************************
  subroutine energy_info()
    !trapezium rule to calculate the integral
    implicit none
    real :: sdot(3), sdoti(3)
    integer :: infront
    integer :: i
    energy=0.
    do i=1, pcount
      if (f(i)%infront==0) cycle !check for 'empty' particles
      infront=f(i)%infront
      call get_deriv_1(i,sdot)
      call get_deriv_1(infront,sdoti)
      energy=energy+0.5*dist_gen(f(i)%x,f(i)%ghosti)*&
      (dot_product(f(i)%u,cross_product(f(i)%x,sdot))+&
       dot_product(f(infront)%u,cross_product(f(infront)%x,sdoti)))
    end do
  end subroutine
  !*************************************************
  subroutine curv_info()
    !caculate the mean, min, max curvature of the vortex system
    implicit none
    real, allocatable :: curvi(:)
    integer :: i, j
    !------------histogram parameters below---------------------
    !warning - at present the matlab routine to plot the histogram 
    !is not adpative therefore if the number of bins is changed 
    !the script must also be changed - warning
    integer, parameter :: bin_num=10 !number of bins
    real :: bin_size, bins(bin_num) !the bins of the histograms
    real :: probability(bin_num) !the probability of that bin
    allocate(curvi(pcount)) !allocate this array pcount size
    do i=1, pcount
      if (f(i)%infront==0) then
        curvi(i)=0. !check for 'empty' particles
      else
        curvi(i)=curvature(i) !general.mod
      end if
    end do
    !compute the average/min/max of this array
    kappa_bar=sum(curvi)/count(mask=f(:)%infront>0)
    kappa_max=maxval(curvi)
    kappa_min=minval(curvi,mask=curvi>0)
    !experimental creation of a histogram
    if (curv_hist) then
      !first we set the bins
      bin_size=(kappa_max-kappa_min)/bin_num
      bins=0.
      !now loop over particles and bins and create the histogram
      do i=1, pcount
        if (f(i)%infront==0) cycle
        do j=1, bin_num
          if (curvi(i)>kappa_min+(j-1)*bin_size) then
            if (curvi(i)<kappa_min+j*bin_size) then
              bins(j)=bins(j)+1
            end if
          end if
        end do
      end do
      bins=bins/(count(mask=f(:)%infront>0)*bin_size) !normalise 
      open(unit=79,file='data/curv_pdf.log',position='append')
      write(79,*) '%----------------t=',t,'---------------'
      do j=1, bin_num
        write(79,*) kappa_min+(j-1)*bin_size, bins(j) 
      end do 
      close(79)
    end if
    deallocate(curvi) !deallocate helper array
  end subroutine
end module
