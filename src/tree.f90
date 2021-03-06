!>all the tree routines used by the code are contained within this module,
!>pointers are used to create a linked list, be very careful to avoid memory
!>leaks!
module tree
   ! CONTAINS ALL THE ROUTINES USED TO EVALUATE THE BIOT-SAVART INTEGRAL
   ! USING TREE CODE METHODS - A WARNING, POINTERS ARE USED HEAVILY HERE
   ! SO BE VERY CAREFUL, MEMORY LEAKS ARE A DISTINCT POSSIBILITY!
   use cdata
   use general
   !>the tree structure, uses pointers
   !>@param posx @param posy @param posz bottem left corner of box
   !>@param width the width of the mesh 
   !>@param fbl child cell front bottom left, @param fbr child cell front bottom right
   !>@param ftl child cell front top left @param ftr child cell front top right
   !>@param bbl child cell back bottom left, @param bbr child cell back bottom right
   !>@param btl child cell back top left @param btr child cell back top right
   !>@param all particles in the cell
   !>@param centx @param centy @param centz  @param cent the center of 'mass' of the cell
   !>@param cent_u center of mass' velocity
   !>@param Cij a sort of circulation tensor
   !>@param circ the total circulation (vector) of the cell
   type node
      real :: posx,posy,posz 
      real :: width 
      integer :: pcount 
      type (node), pointer :: fbl, fbr !front (y) bottom left, right (child)
      type (node), pointer :: ftl, ftr !front (y) top left, right (child)
      type (node), pointer :: bbl, bbr !back (y) bottom left, right (child)
      type (node), pointer :: btl, btr !back (y) top left, right (child)
      type (qvort), allocatable :: parray(:) 
      real :: centx, centy, centz, cent(3)
      real :: cent_u(3) !velocity of the center (\sum u_j)
      real :: circ(3) !total circulation vector
      real :: Cij(3, 3) ! \int r_i dr_j where r are relative to the center of mass
                        ! a "circulation tensor", for a lack of better name
   end type node
   !>the vortex point tree
   type (node), pointer :: vtree 
   !>how many evaluations are required to calculate the tree vel field
   integer, protected :: eval_counter
  contains
  !>An array to set the main node of the tree and then call the recurisive 
  !>subroutine build_tree this is called in the main program
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
    vtree%cent = [vtree%centx, vtree%centy, vtree%centz]
    vtree%circ=0. !0 the circulation
    vtree%Cij = 0.
    !now nullify its children - for safety
    nullify(vtree%fbl, vtree%fbr, vtree%ftl, vtree%ftr)
    nullify(vtree%bbl, vtree%bbr, vtree%btl, vtree%btr)
    call build_tree(vtree) !the recursive routine which builds the tree
    !set the evaluation counter here
    select case(velocity)
      case('LIA')
        eval_counter=count(mask=f(:)%infront>0)
      case('BS')
        if (periodic_bc) then
          eval_counter=27*count(mask=f(:)%infront>0)**2
        else if (periodic_bc_notx) then
          eval_counter=9*count(mask=f(:)%infront>0)**2
        else if (periodic_bc_notxy) then
          eval_counter=3*count(mask=f(:)%infront>0)**2
        else
          eval_counter=count(mask=f(:)%infront>0)**2
        end if
      case('Tree')
        eval_counter=2*count(mask=f(:)%infront>0) !set to 2 as we ignore particle/particle behind
    end select
  end subroutine
  !**************************************************************************
  !>called by above, build the tree, continue to subdivide until a cell is empty
  !>or contains only one particle - octree structure used each new cell is created
  !>by the routine below new_tree
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
      if ((mod(itime,shots)==0).and.tree_print) then
        write(unit=print_file,fmt="(a,i3.3,a)")"./data/tree_mesh",itime/shots,".log"
        open(unit=97,file=print_file,action='write', position='append')
          write(97,*) vtree%posx, vtree%posx+vtree%width, vtree%posy, &
                     vtree%posy+vtree%width, vtree%posz, vtree%posz+vtree%width, vtree%pcount
        close(97)
      end if
    end if
  end subroutine
  !********************************************************************************
  !>allocate a new cell, populate with particles and calculate centre of mass, 
  !>total circulation vector
  subroutine new_tree(vtree,posx,posy,posz,width, loc_pcount,parray)
    !create a new tree
    implicit none
    integer, intent(IN) :: loc_pcount
    real, intent(IN) :: posx, posy, posz, width
    type(qvort) :: parray(loc_pcount)
    type(node), pointer :: vtree
    integer :: i, counter, k
    allocate(vtree)
    allocate(vtree%parray(loc_pcount))
    !set up this node of the tree with the input info
    vtree%posx=posx ; vtree%posy=posy ; vtree%posz=posz
    vtree%width=width ; vtree%centx=0. ; vtree%centy=0.
    vtree%centz=0. ; vtree%circ=0. ; vtree%cent_u = 0.
    vtree%Cij = 0
    vtree%parray(:)%infront=0 !0 the infront marker
    counter=0 !0 the counter
    do i=1, loc_pcount
      if (parray(i)%infront==0) cycle
      if((parray(i)%x(1)>=posx).and.(parray(i)%x(1)<(posx+width)).and. &
         (parray(i)%x(2)>=posy).and.(parray(i)%x(2)<(posy+width)).and. &
         (parray(i)%x(3)>=posz).and.(parray(i)%x(3)<(posz+width))) then
        counter=counter+1
        vtree%parray(counter)=parray(i)
        !add the position of the particle to the position of the mesh
        vtree%centx=vtree%centx+parray(i)%x(1)
        vtree%centy=vtree%centy+parray(i)%x(2)
        vtree%centz=vtree%centz+parray(i)%x(3)
        !centre's velocity is the mean of velocities
        vtree%cent_u = vtree%cent_u + parray(i)%u
        !also account for the circulation
        vtree%circ = vtree%circ + parray(i)%ghosti - parray(i)%x
        do k=1, 3
           vtree%Cij(:,k) = vtree%Cij(:,k)&
                &         + parray(i)%x*(parray(i)%ghosti(k) - parray(i)%x(k))
        end do
      end if
    end do
    !take the average of the positions
    if (counter/=0) then
     vtree%centx=vtree%centx/counter
     vtree%centy=vtree%centy/counter
     vtree%centz=vtree%centz/counter
     vtree%cent_u = vtree%cent_u / counter

     vtree%cent = [vtree%centx, vtree%centy, vtree%centz] !it's easier to work with vectors
     do i=1, 3
        do k=1, 3
           vtree%Cij(i, k) = vtree%Cij(i, k) - vtree%cent(i)*vtree%circ(k)
        end do
     end do
    end if
    vtree%pcount=counter
    !nullify the children
    nullify(vtree%fbl, vtree%fbr, vtree%ftl, vtree%ftr)
    nullify(vtree%bbl, vtree%bbr, vtree%btl, vtree%btr)
    call build_tree(vtree) !recall the routine that builds the tree
  end subroutine
  !********************************************************************************
  !>empty out the tree structure to avoid a memory leak - repeatedly calls dead_tree
  recursive subroutine empty_tree(vtree)
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
  !>deallocate a tree
  subroutine dead_tree(vtree)
    !chop down that tree (metaphorically) 
    implicit none
    type (node), pointer :: vtree
    call empty_tree(vtree) !kill its children (if it has any)
    !print*, vtree%pcount
    deallocate(vtree%parray)
    deallocate(vtree)
    nullify(vtree)
  end subroutine
  !********************************************************************************
  !>find the closest particle using the tree structure, theoretically should be NlogN
  !>repeatedly calls closest_tree
  subroutine pclose_tree
    !loop over all the particles and find the closest one to i using the tree mesh
    implicit none
    integer :: i
    !$omp parallel do private(i)
    do i=1, pcount
      f(i)%closestd=100. !arbitrarily high
      if (f(i)%infront==0) cycle !empty particle
      if (mirror_bc) then
        !do not test if you are pinned 
        if (f(i)%pinnedi.or.f(i)%pinnedb) cycle 
      end if 
      call closest_tree(i,vtree)
    end do
    !$omp end parallel do
  end subroutine
  !********************************************************************************
  !>recurisve routine used by pclose_tree to find the closest particle to i
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
           if (mirror_bc) then
             !test to see if j is pinned 
             if (f(j)%pinnedi.or.f(j)%pinnedb) then
               !do nothing its pinned
             else
               f(i)%closest=j ; f(i)%closestd=dist !OK not pinned
             end if 
           else
             f(i)%closest=j ; f(i)%closestd=dist
           end if
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
  !********************************************************************************
  !>find the closest particle using the tree structure
  !>this particle cannot be on the same loop, theoretically should be NlogN
  !>repeatedly calls closest_tree, used by line_sep diagnostic routine
  subroutine pclose_tree_loop
    !loop over all the particles and find the closest one to i using the tree mesh
    implicit none
    integer :: i
    !$omp parallel do private(i)
    do i=1, pcount
      f(i)%closestd_loop=100. !arbitrarily high
      if (f(i)%infront==0) cycle !empty particle
      call closest_tree_loop(i,vtree)
    end do
    !$omp end parallel do
  end subroutine
  !********************************************************************************
  !>recurisve routine used by pclose_tree_loop to find the closest particle to i
  recursive subroutine closest_tree_loop(i,vtree)
    !find the nearest particle to i who is not on the same loop
    implicit none
    integer, intent(IN) :: i
    type (node), pointer :: vtree
    real :: theta, dist !tree distance
    integer :: j
    logical :: same_loop
     if (f(i)%infront==0) return !zero particle exit
     if (vtree%pcount==0) return !empty box no use
     !work out distances opening angles etc. here
     dist=sqrt((f(i)%x(1)-vtree%centx)**2+(f(i)%x(2)-vtree%centy)**2+(f(i)%x(3)-vtree%centz)**2)
     theta=vtree%width/dist
     if (vtree%pcount==1.or.theta<.1) then !fairly tight critical opening angle
       !see if the distance if less than the minimum distance
       if (dist<f(i)%closestd_loop) then
         j=f(vtree%parray(1)%infront)%behind
         if ((vtree%pcount==1).and.(i/=j).and.(f(i)%infront/=j).and.(f(i)%behind/=j)) then
           !the above line ignores contribution from neighbouring particles
           !ensuring there is only one particle in the box
           call same_loop_test(i,j,same_loop)
           if (same_loop.eqv..false.) then
             f(i)%closestd_loop=dist
           end if
         end if
       end if
     else
       !open the box up and use the child cells ; make sure its the right routine!
       call closest_tree_loop(i,vtree%fbl) ; call closest_tree_loop(i,vtree%fbr)
       call closest_tree_loop(i,vtree%ftl) ; call closest_tree_loop(i,vtree%ftr)
       call closest_tree_loop(i,vtree%bbl) ; call closest_tree_loop(i,vtree%bbr)
       call closest_tree_loop(i,vtree%btl) ; call closest_tree_loop(i,vtree%btr)
     end if
   end subroutine
   !************************************************************************
   !>get the tree code approximation to the biot savart integral, giving
   !>the induced velocity at particle i due to all other vortices, accepts 
   !>an arguement shift, which allows you to shift the position of the
   !>vortices for periodic boundary conditions
   recursive subroutine tree_walk(i,vtree,shift,u)
     implicit none
     integer, intent(IN) :: i !the particle we are interested in
     real, intent(IN) :: shift(3) !shift tree (periodicity)
     type (node), pointer :: vtree !the tree
     real, intent(inout) :: u(3) !the velocity at particle i
     real :: vect(3) !helper vector
     real :: dist, theta !helper variables
     real :: geo_dist !distance between centre of vorticity and centre of cell
     real :: a_bs, b_bs, c_bs, u_bs(3) !helper variables 
     integer :: j=0 !the particle in the tree mesh
     integer :: k, l
     if(f(i)%infront == 0) return
     if (vtree%pcount==0) return !empty box no use
     !work out distances opening angles etc. here
     dist=sqrt((f(i)%x(1)-(vtree%centx+shift(1)))**2+&
         &     (f(i)%x(2)-(vtree%centy+shift(2)))**2+&
         &     (f(i)%x(3)-(vtree%centz+shift(3)))**2)
     if (tree_extra_correction) then
        geo_dist=sqrt((vtree%centx-(vtree%posx+vtree%width/2.))**2+&
                &     (vtree%centy-(vtree%posy+vtree%width/2.))**2+&
                &     (vtree%centz-(vtree%posz+vtree%width/2.))**2)
        if (dist-geo_dist>0.) then
           theta=vtree%width/(dist-geo_dist)
        else
           theta=vtree%width/dist 
        end if
     else
        theta=vtree%width/dist 
     end if
     if (vtree%pcount==1) then
        !check that the particle is not itself or the particle behind
        j=f(vtree%parray(1)%infront)%behind
        if (j==i) return
        if (j==f(i)%behind) return
        if (j==f(i)%infront) return
        if (dist<epsilon(0.)) then
           call fatal_error('tree.mod:tree_walk', & 
                'singularity in BS (tree) velocity field - &
                this is normally caused by having recon_shots too large') !cdata.mod
        end if
        !add the contribution of single segment
        !vtree%cent is the same as vtree%parray(1)%x and
        !vtree%circ is == f(vtree%parray(1)%infront)%x
        a_bs = vector_norm(vtree%parray(1)%x + shift - f(i)%x) !lRj
        b_bs = vector_norm(vtree%parray(1)%ghosti + shift - f(i)%x) !lRjp1

        !denom
        c_bs = a_bs*b_bs*(a_bs*b_bs&
             &           + dot_product(vtree%parray(1)%x + shift - f(i)%x,&
             &                         vtree%parray(1)%ghosti + shift - f(i)%x))

        u_bs = (a_bs + b_bs)/c_bs*cross_product(vtree%parray(1)%x + shift - f(i)%x,&
             &                                  vtree%parray(1)%ghosti + shift - f(i)%x)
        u_bs = u_bs * quant_circ/4/pi

        u = u + u_bs
        eval_counter=eval_counter+1 !increment this counter
     else if (theta < tree_theta) then !use the contribution of this cell
        !there is more than one vortex particle in this node so use
        !the first two terms in Taylor expansion of the Biot-Savart
        !integral
        !vect = vtree%cent + shift - f(i)%x

        vect = f(i)%x - (vtree%cent + shift)

        u_bs = cross_product(vtree%circ, vect)/dist**3;
        do k=1, 3
           do l=1, 3
              u_bs = u_bs - eps(:, k, l)*vtree%Cij(k, l)/dist**3;
              u_bs = u_bs + 3*eps(:, k, l)*vect(k)*dot_product(vect, vtree%Cij(:,l))&
                   &        /dist**5
           end do
        end do

        u = u + u_bs*quant_circ/4./pi

        eval_counter=eval_counter+1
     else
        !open the box up and use the child cells
        call tree_walk(i,vtree%fbl,shift,u) ; call tree_walk(i,vtree%fbr,shift,u)
        call tree_walk(i,vtree%ftl,shift,u) ; call tree_walk(i,vtree%ftr,shift,u)
        call tree_walk(i,vtree%bbl,shift,u) ; call tree_walk(i,vtree%bbr,shift,u)
        call tree_walk(i,vtree%btl,shift,u) ; call tree_walk(i,vtree%btr,shift,u)
     end if
   end subroutine
!************************************************************************
   !>a more general form of above, where a position is an input
   !>this is used to calculate the superfluid velocity at mesh points
   recursive subroutine tree_walk_general(x,vtree,shift,u)
     implicit none
     real, intent(IN) :: x(3) !the position we want the velocity at
     real, intent(IN) :: shift(3) !shift tree (periodicity)
     type (node), pointer :: vtree !the tree
     real :: u(3) !the velocity at particle i
     real :: vect(3) !helper vector
     real :: a_bs, b_bs, c_bs, u_bs(3) !helper variables 
     real :: dist, theta !helper variables
     integer :: j=0 !the particle in the tree mesh
     integer :: k,l

     if (vtree%pcount==0) return !empty box no use
     !work out distances opening angles etc. here
     dist=sqrt((x(1)-(vtree%centx+shift(1)))**2+&
          (x(2)-(vtree%centy+shift(2)))**2+&
          (x(3)-(vtree%centz+shift(3)))**2)
     if (dist<epsilon(0.)) return !avoid 1/0.

     theta=vtree%width/dist !the most simple way we can improve this
     if (vtree%pcount==1) then
        !add the contribution of single segment
        !vtree%cent is the same as vtree%parray(1)%x and
        !vtree%circ is == f(vtree%parray(1)%infront)%x
        a_bs = vector_norm(vtree%cent + shift - x) !lRj
        b_bs = vector_norm(vtree%cent + vtree%circ + shift - x) !lRjp1

        !denom
        c_bs = a_bs*b_bs*(a_bs*b_bs&
             &           + dot_product(vtree%cent + shift - x,&
             &                         vtree%cent + vtree%circ + shift - x))

        u_bs = (a_bs + b_bs)/c_bs*cross_product(vtree%cent + shift - x,&
             &                                  vtree%cent + vtree%circ + shift - x)

        u = u + u_bs/4/pi
        eval_counter=eval_counter+1 !increment this counter

     else if (theta < tree_theta) then !use the contribution of this cell
        vect=x-(vtree%cent+shift)
        u_bs = cross_product(vtree%circ, vect)/dist**3;
        do k=1, 3
           do l=1, 3
              u_bs = u_bs - eps(:, k, l)*vtree%Cij(k, l)/dist**3;
              u_bs = u_bs + 3*eps(:, k, l)*vect(k)*dot_product(vect, vtree%Cij(:,l))&
                   &        /dist**5
           end do
        end do

        if (hollow_mesh_core) then
           if (dist<delta) then
              u_bs=0. !0 the velocity
           end if
        end if
        u=u+u_bs
        eval_counter=eval_counter+1 !increment this counter
     else
        !open the box up and use the child cells
        call tree_walk_general(x,vtree%fbl,shift,u) ; call tree_walk_general(x,vtree%fbr,shift,u)
        call tree_walk_general(x,vtree%ftl,shift,u) ; call tree_walk_general(x,vtree%ftr,shift,u)
        call tree_walk_general(x,vtree%bbl,shift,u) ; call tree_walk_general(x,vtree%bbr,shift,u)
        call tree_walk_general(x,vtree%btl,shift,u) ; call tree_walk_general(x,vtree%btr,shift,u)
     end if
   end subroutine

   subroutine tree_mat_diff(vtree, x, dtvs, vs_g_vs, vs, periodic, cutoff)
     type(node), intent(in)        :: vtree
     real, intent(in)              :: x(3)
     real, intent(out)             :: dtvs(3), vs_g_vs(3), vs(3)
     real, intent(in), optional    :: cutoff
     logical, intent(in), optional :: periodic

     real :: vec(3), dist, theta
     !periodicity
     integer :: px_l = 0, py_l = 0, pz_l = 0
     integer :: px_h = 0, py_h = 0, pz_h = 0
     integer :: px, py, pz
     real    :: shift(3)

     if(present(periodic) .and. periodic) then
        px_l = -1; py_l = -1; pz_l = -1
        px_h =  1; py_h =  1; pz_h =  1
     end if

     if(vtree%pcount == 0) return;
     do px=px_l, px_h ; do py=py_l, py_h; do pz=pz_l, pz_h
        shift = [px, py, pz]*box_size

     end do; end do; end do
     
   end subroutine tree_mat_diff
end module
