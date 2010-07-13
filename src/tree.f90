module tree
   ! CONTAINS ALL THE ROUTINES USED TO EVALUATE THE BIOT-SAVART INTEGRAL
   ! USING TREE CODE METHODS - A WARNING, POINTERS ARE USED HEAVILY HERE
   ! SO BE VERY CAREFUL, MEMORY LEAKS ARE A DISTINCT POSSIBILITY!
   use cdata
   use general
   type node
      real :: posx,posy,posz !bottem left corner of box
      real :: width !width (obviously)
      integer :: pcount !particle count in this box
      type (node), pointer :: fbl, fbr !front (y) bottom left, right (child)
      type (node), pointer :: ftl, ftr !front (y) top left, right (child)
      type (node), pointer :: bbl, bbr !back (y) bottom left, right (child)
      type (node), pointer :: btl, btr !back (y) top left, right (child)
      type (qvort), allocatable :: parray(:) !the particles in that mesh
      real :: centx, centy, centz !centres of 'mass'
      real :: circ(3) !total circulation vector
   end type node
   type (node), pointer :: vtree  ! our vortex tree
  contains
  subroutine construct_tree()
    allocate(vtree)
    allocate(vtree%parray(pcount))
    !set up the main node, not a 100% sure if this needs doing as will never 
    !be used for evaluations
    vtree%posx=-box_size/2. ; vtree%posy=-box_size/2. ; vtree%posz=-box_size/2.
    vtree%width=box_size
    vtree%pcount=pcount
    vtree%parray=f
    vtree%centx=sum(f(:)%x(1), mask=f(:)%infront>0)/count(mask=f(:)%infront>0)
    vtree%centy=sum(f(:)%x(2), mask=f(:)%infront>0)/count(mask=f(:)%infront>0)
    vtree%centz=sum(f(:)%x(3), mask=f(:)%infront>0)/count(mask=f(:)%infront>0)
    !now nullify its children - for safety
    nullify(vtree%fbl, vtree%fbr, vtree%ftl, vtree%ftr)
    nullify(vtree%bbl, vtree%bbr, vtree%btl, vtree%btr)
    call build_tree(vtree) !the recursive routine which builds the tree
  end subroutine
  !**************************************************************************
  recursive subroutine build_tree(vtree)
    implicit none
    type(node), pointer :: vtree
    character (len=40) :: print_file
    if (vtree%pcount>1) then
      !divide up into 8 children
      !front, bottom, left (x,y,z)
      call new_tree(vtree%fbl, vtree%posx, vtree%posy, vtree%posz, &
                    vtree%width/2., vtree%pcount, vtree%parray)
      !front, bottom, right (x+w,y,z)
      call new_tree(vtree%fbr, vtree%posx+(vtree%width/2), vtree%posy, vtree%posz, &
                    vtree%width/2., vtree%pcount, vtree%parray)
      !front, top, left (x,y,z+w)
      call new_tree(vtree%ftl, vtree%posx, vtree%posy,vtree%posz+(vtree%width/2), &
                    vtree%width/2., vtree%pcount, vtree%parray)
      !front, top, right (x+w,y,z+w)
      call new_tree(vtree%ftr, vtree%posx+(vtree%width/2), vtree%posy,vtree%posz+(vtree%width/2), &
                    vtree%width/2., vtree%pcount, vtree%parray)
      !back, bottom, left (x,y+w,z)
      call new_tree(vtree%bbl, vtree%posx, vtree%posy+(vtree%width/2), vtree%posz, &
                    vtree%width/2., vtree%pcount, vtree%parray)
      !back, bottom, right (x+w,y+w,z)
      call new_tree(vtree%bbr, vtree%posx+(vtree%width/2), vtree%posy+(vtree%width/2), vtree%posz, &
                    vtree%width/2., vtree%pcount, vtree%parray)
      !back, top, left (x,y+w,z+w)
      call new_tree(vtree%btl, vtree%posx, vtree%posy+(vtree%width/2), vtree%posz+(vtree%width/2), &
                    vtree%width/2., vtree%pcount, vtree%parray)
      !back, top, right (x+w,y+w,z+w)
      call new_tree(vtree%btr, vtree%posx+(vtree%width/2), vtree%posy+(vtree%width/2), vtree%posz+(vtree%width/2), &
                    vtree%width/2., vtree%pcount, vtree%parray)
    else
      !print to file to test the algorithm works
      if (mod(itime,shots)==0) then
        write(unit=print_file,fmt="(a,i3.3,a)")"./data/tree_mesh",itime/shots,".log"
        open(unit=97,file=print_file,action='write', position='append')
          write(97,*) vtree%posx, vtree%posx+vtree%width, vtree%posy, &
                     vtree%posy+vtree%width, vtree%posz, vtree%posz+vtree%width, vtree%pcount
        close(97)
      end if
    end if
  end subroutine
  !********************************************************************************
  subroutine new_tree(vtree,posx,posy,posz,width, loc_pcount,parray)
    !create a new tree
    implicit none
    integer, intent(IN) :: loc_pcount
    real, intent(IN) :: posx, posy, posz, width
    type(qvort) :: parray(loc_pcount)
    type(node), pointer :: vtree
    integer :: i, counter
    allocate(vtree)
    allocate(vtree%parray(loc_pcount))
    !set up this node of the tree with the input info
    vtree%posx=posx ; vtree%posy=posy ; vtree%posz=posz
    vtree%width=width ; vtree%centx=0. ; vtree%centy=0.
    vtree%centz=0. ; vtree%circ=0.
    vtree%parray(:)%infront=0 !0 the infront marker
    counter=0 !0 the counter
    do i=1, loc_pcount
      if (parray(i)%infront==0) cycle
      if((parray(i)%x(1)>=posx).and.(parray(i)%x(1)<(posx+width)).and. &
         (parray(i)%x(2)>=posy).and.(parray(i)%x(2)<(posy+width)).and. &
         (parray(i)%x(3)>=posz).and.(parray(i)%x(3)<(posz+width))) then
        counter=counter+1
        vtree%parray(counter)=parray(i)
        vtree%centx=vtree%centx+parray(i)%x(1)
        vtree%centy=vtree%centy+parray(i)%x(2)
        vtree%centz=vtree%centz+parray(i)%x(3)
      end if
    end do
    !take the average of the positions
    if (counter/=0) then
     vtree%centx=vtree%centx/counter
     vtree%centy=vtree%centy/counter
     vtree%centz=vtree%centz/counter
    end if
    vtree%pcount=counter
    !nullify the children
    nullify(vtree%fbl, vtree%fbr, vtree%ftl, vtree%ftr)
    nullify(vtree%bbl, vtree%bbr, vtree%btl, vtree%btr)
    call build_tree(vtree) !recall the routine that builds the tree
  end subroutine
  !********************************************************************************
  recursive subroutine empty_tree(vtree)
    !empty out the tree structure to avoid a memory leak
    implicit none
    type (node), pointer :: vtree
    if(vtree%pcount>1) then
      !kill all it's 8 children
      call dead_tree(vtree%fbl)
      call dead_tree(vtree%fbr)
      call dead_tree(vtree%ftl)
      call dead_tree(vtree%ftr)
      call dead_tree(vtree%bbl)
      call dead_tree(vtree%bbr)
      call dead_tree(vtree%btl)
      call dead_tree(vtree%btr)
    end if
  end subroutine
  !********************************************************************************
  subroutine dead_tree(vtree)
    !chop down that tree (metaphorically) 
    implicit none
    type (node), pointer :: vtree
    call empty_tree(vtree) !kill its children (if it has any)
    deallocate(vtree%parray)
    deallocate(vtree)
    nullify(vtree)
  end subroutine
  !********************************************************************************
  subroutine pclose_tree
    !loop over all the particles and find the closest one to i using the tree mesh
    implicit none
    integer :: i
    do i=1, pcount
      f(i)%closestd=10. !arbitrarily high
      call closest_tree(i,vtree)
    end do
  end subroutine
  !********************************************************************************
  recursive subroutine closest_tree(i,vtree)
    !find the nearest particle to i
    implicit none
    integer, intent(IN) :: i
    type (node), pointer :: vtree
    real :: theta, dist !tree distance
    integer :: j
     if (f(i)%infront==0) return !zero particle exit
     if (vtree%pcount==0) return !empty box no use
     !work out distances opening angles etc. here
     dist=sqrt((f(i)%x(1)-vtree%centx)**2+(f(i)%x(2)-vtree%centy)**2+(f(i)%x(3)-vtree%centz)**2)
     theta=vtree%width/dist
     if (vtree%pcount==1.or.theta<.5) then
       !see if the distance if less than the minimum distance
       if (dist<f(i)%closestd) then
         j=f(vtree%parray(1)%infront)%behind
         if (vtree%pcount>1) then
           !skip this the particle must be isolated
           !probably will not happen in practice
           f(i)%closest=0
         else if ((i/=j).and.(f(i)%infront/=j).and.(f(i)%behind/=j)) then
           !the above line ensures we do not reconnect with particles infront/behind
           f(i)%closest=j ; f(i)%closestd=dist
         end if
       end if
     else
       !open the box up and use the child cells
       call closest_tree(i,vtree%fbl) ; call closest_tree(i,vtree%fbr)
       call closest_tree(i,vtree%ftl) ; call closest_tree(i,vtree%ftr)
       call closest_tree(i,vtree%bbl) ; call closest_tree(i,vtree%bbr)
       call closest_tree(i,vtree%btl) ; call closest_tree(i,vtree%btr)
     end if
   end subroutine
end module
