c- constants -----------------------------------------------------------

! #define tSTATSTART uparam(1) /* start time for averaging */
! #define tSTATFREQ  uparam(2) /* output frequency for statistics */

c data extraction along wall normal direction
! #define INTP_NMAX 200 /* number of sample points */
! #define XCINT 1.0     /* x coordinate of 1D line*/
! #define ZCINT 1.0     /* z coordinate of 1D line */

c mesh dimensions
! #define BETAM 2.4     /* wall normal stretching parameter */
! #define PI (4.*atan(1.))
! #define XLEN (2.*PI)
! #define ZLEN PI
! #define NUMBER_ELEMENTS_X 16
! #define NUMBER_ELEMENTS_Y 12
! #define NUMBER_ELEMENTS_Z 8

c-----------------------------------------------------------------------
      subroutine uservp (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      integer e

      utrans = 1.
      udiff  = param(2)

      if (ifield .eq. 2) then
         e = gllel(ieg)
         udiff = param(8)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine userf  (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      ffx = 0.01
      ffy = 0.0
      ffz = 0.0

      return
      end
c-----------------------------------------------------------------------
      subroutine userq  (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      qvol =  0.0

      return
      end
c-----------------------------------------------------------------------
      subroutine userchk

      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'SOLN'
      include 'MVGEOM'
      include 'INPUT'
      include 'TSTEP'
!      include 'TOTAL'
      include 'MASS'

      include 'F3D'
      include 'FS_ALE'

      include 'TEST'

      integer ntot1,ntot2
      integer i,j

      integer igeom
      character cb*3
      integer ie,iface,nfaces

      if (istep.eq.0) then
        call frame_start
      endif  

      call frame_monitor

!      call slp_mark_faces()

      ifto = .true.

      ntot1 = lx1*ly1*lz1*nelv
      ntot2 = lx2*ly2*lz2*lelv


      if (istep.eq.0) then

        call initp_f3d
        call rzero(tmp4,ntot2)

      endif

      call copy(t,vz,ntot1)

      call test_random

      call fs_mvmesh()

      if (istep.eq.nsteps.or.lastep.eq.1) then
        call frame_end
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine userbc (ix,iy,iz,iside,ieg)

      implicit none

      include 'SIZE'
!      include 'TOTAL'
      include 'NEKUSE'
      include 'GEOM'

      integer ix,iy,iz,iside,ieg
      real pi

      integer jp
      common /ppointr/ jp

      real r1,r2,omega1,omega2
      real a1,a2
      save r1,r2,a1,a2,omega1,omega2
      real rmid

      integer isave
      save isave
      data isave /0/

      real glmin,glmax
      integer n

      
      if (isave.eq.0) then
        n = lx1*ly1*lz1*nelv
        r1 = glmin(ym1,n)
        r2 = glmax(ym1,n)
        omega1 = 1.0/r1
        omega2 = 0.0
        a1 = (omega2*(r2**2) - omega1*r1*r1)/(r2**2 - r1**2)
        a2 = (omega1 - omega2)*(r1**2)*(r2**2)/(r2**2 - r1**2)
        isave = 1
      endif

      if (jp.eq.0) then
        ux   = 0.0
        uy   = 0.0
        uz   = 0.0
        
        rmid = (r1+r2)/2
        if (y.lt.(rmid)) uz   = 1.0
!        temp = 0.0
!        temp = a1*y + a2/y
      else
        ux   = 0.0
        uy   = 0.0
        uz   = 0.0
        temp = 0.0
      endif  

      return
      end
c-----------------------------------------------------------------------
      subroutine useric (ix,iy,iz,ieg)

      implicit none

      include 'SIZE'
      include 'INPUT'         ! if3d
      include 'PARALLEL'
      include 'NEKUSE'
      include 'GEOM'

      include 'F3D'

      integer ix,iy,iz,ieg
      real pi

      integer jp
      common /ppointr/ jp

      real fcoeff(3)
      real xl(3)
      real mth_ran_dst

      logical ifcouette
      logical ifpoiseuille
      logical iftaylor
      logical iftestmvb

      real r1,r2,omega1,omega2
      real a1,a2
      save r1,r2,a1,a2,omega1,omega2

      integer isave
      save isave
      data isave /0/

      real glmin,glmax
      integer n


      ifcouette         = .false.
      ifpoiseuille      = .false.
      iftaylor          = .false.
      iftestmvb         = .true.

      pi = 4.0*atan(1.0)

      if (isave.eq.0) then
        n = lx1*ly1*lz1*nelv
        r1 = glmin(ym1,n)
        r2 = glmax(ym1,n)
        omega1 = 1.0/r1
        omega2 = 0.0
        a1 = (omega2*(r2**2) - omega1*r1*r1)/(r2**2 - r1**2)
        a2 = (omega1 - omega2)*(r1**2)*(r2**2)/(r2**2 - r1**2)
        isave = 1
      endif


      if (jp.eq.0) then

        if (ifpoiseuille) then
          ux = 1.0 - y**2
          uy = 0.
          uz = 0.0 + 0.0
        elseif (ifcouette) then
          ux = 0.0 + 1.0*y
          uy = 0.
          uz = 0.0 + 0.0
        elseif (iftaylor) then
!          ux = 0.0 + 0.0
!          uy = 0.
!          uz = a1*y + a2/y
!         Testing  
          temp = 0.01*(a2*y + a1/y)   ! just testing
          ux   = temp
          uy   = temp
        elseif (iftestmvb) then
          ux   = -0.05*exp(-((y-r1)/0.25)**2)
          uy   = -0.01*exp(-((y-r1)/0.25)**2)
          uz   = 0.1*(a1*y + a2/y)   ! just testing
!          temp = 0.01*(a1*y + a2/y)   ! just testing
        endif  
      else

!       perturbation; white noise
        xl(1) = X
        xl(2) = Y
        if (IF3D.or.iff3d) xl(3) = Y+X
        
        if (jp.eq.1) then
          fcoeff(1)=  3.0e4
          fcoeff(2)= -1.5e3
          fcoeff(3)=  0.5e5
        else
          fcoeff(1)=  9.0e4
          fcoeff(2)=  1.5e3
          fcoeff(3)= -2.5e5
        endif          
        ux=UPARAM(1)*mth_ran_dst(ix,iy,iz,ieg,xl,fcoeff)
        if (jp.eq.1) then
          fcoeff(1)=  2.3e4
          fcoeff(2)=  2.3e3
          fcoeff(3)= -2.0e5
        else
          fcoeff(1)=  1.3e4
          fcoeff(2)= -5.8e3
          fcoeff(3)= -1.9e5
        endif
        uy=UPARAM(1)*mth_ran_dst(ix,iy,iz,ieg,xl,fcoeff)
        if (IF3D.or.iff3d) then
           if (jp.eq.1) then           
             fcoeff(1)= 2.e4
             fcoeff(2)= 1.e3
             fcoeff(3)= 1.e5
           else
             fcoeff(1)= -1.9e4
             fcoeff(2)= -8.0e3
             fcoeff(3)=  3.2e5
           endif 
           uz=UPARAM(1)*mth_ran_dst(ix,iy,iz,ieg,xl,fcoeff)
        else
           uz = 0.0
        endif

      endif


      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat   ! This routine to modify element vertices

      implicit none  
  
      include 'SIZE'      ! _before_ mesh is generated, which 
      include 'INPUT'
      include 'GEOM'
      include 'TSTEP'
!      include 'TOTAL'     ! guarantees GLL mapping of mesh.

      integer n,i,j
      real r0

!      ifaxis = .true.   ! just for initialization
      param(42)=1       ! 0: GMRES (nonsymmetric), 1: PCG w/o weights
      param(43)=1       ! 0: Additive multilevel (param 42=0), 1: Original 2 level
      param(44)=1       ! 0: E based Schwartz, 1: A based Schwartz

      n = nelv * 2**ldim
!      xmin = glmin(xc,n)
!      xmax = glmax(xc,n)
!      ymin = glmin(yc,n)
!      ymax = glmax(yc,n)
!      zmin = glmin(zc,n)
!      zmax = glmax(zc,n)
!
!      xscale = XLEN/(xmax-xmin)
!      yscale = 1./(ymax-ymin)
!      zscale = ZLEN/(zmax-zmin)

      pi = 4.*atan(1.0)

      if (abs(uparam(3)).gt.1.0e-6) then
        r0   = abs(uparam(3))
      endif

      if (nio.eq.0) write(6,*) 'R0:', r0

      do j=1,nelv
      do i=1,2**ldim
!         xc(i,j) = 2.0*pi*(xc(i,j))*alphai
         yc(i,j) = yc(i,j) + r0
!         yc(i,1) = tanh(BETAM*(2*yc(i,1)-1))/tanh(BETAM)
!         zc(i,1) = zscale*zc(i,1)
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat2   ! This routine to modify mesh coordinates


      include 'SIZE'
      include 'TOTAL'


!      call outpost(vx,vy,vz,pr,t,'   ')

!      do iel=1,nelt
!      do ifc=1,2*ndim
!         if (cbc(ifc,iel,1) .eq. 'W  ') boundaryID(ifc,iel) = 1 
!         cbc(ifc,iel,2) = cbc(ifc,iel,1) 
!         if (cbc(ifc,iel,1) .eq. 'W  ') cbc(ifc,iel,2) = 't  '
!      enddo
!      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat3
      include 'SIZE'
      include 'TOTAL'

!      param(54) = -1  ! use >0 for const flowrate or <0 bulk vel
                      ! flow direction is given by (1=x, 2=y, 3=z) 
!      param(55) = 1.0 ! flowrate/bulk-velocity 

      call gen_mapping_mvb

      return
      end
c-----------------------------------------------------------------------
!=======================================================================
!> @brief Register user specified modules
      subroutine frame_usr_register
      implicit none

      include 'SIZE'
      include 'FRAMELP'

!     register modules
      call io_register
      call chkpt_register
      call frame_register_f3d
      call tst_register

      return
      end subroutine
!======================================================================
!> @brief Initialise user specified modules
      subroutine frame_usr_init
      implicit none

      include 'SIZE'
      include 'FRAMELP'

!     initialise modules
      call chkpt_init
      call frame_get_param_f3d
      call tst_init

      return
      end subroutine
!======================================================================
!> @brief Finalise user specified modules
      subroutine frame_usr_end
      implicit none

      include 'SIZE'
      include 'FRAMELP'

      
      return
      end subroutine

!-----------------------------------------------------------------------

      subroutine check_vbasea

      implicit none

      include 'SIZE'
      include 'ARN_ARPD'
      include 'TSTEPPERD'
      include 'F3D'
      include 'TSTEP'
      include 'SOLN'
      include 'INPUT'

      integer i,j

      do j=1,2
        i = 1
        call copy(vxp(1,2),vbasea(i,j),tst_nv)
        i = i + tst_nv
        call copy(vyp(1,2),vbasea(i,j),tst_nv)
        i = i + tst_nv
        if (if3d.or.iff3d) then
          call copy(vzp(1,2),vbasea(i,j),tst_nv)
          i = i + tst_nv
        endif
        if (arna_ifpr) then
          call copy(prp(1,2),vbasea(i,j),tst_np)
          i = i + tst_np
        endif  

        ifto = .true.
        call outpost(vxp(1,2),vyp(1,2),vzp(1,2),prp(1,2),
     $               vzp(1,2),'vba')
      enddo 

      call exitt

      return
      endsubroutine check_vbasea
!---------------------------------------------------------------------- 

      subroutine test_random

      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'SOLN'
      include 'INPUT'
      include 'TSTEP'
      include 'MASS'

      include 'F3D'

      include 'TEST'
      include 'IXYZ'
      include 'DXYZ'
      include 'WZ'
      include 'PARALLEL'

      integer ntot1,ntot2
      integer i,j

      integer igeom
      character cb*3
      integer ie,iface,nfaces

      integer nface,nedge,ncrnr
      integer ntotf,ntots,ntotc

      integer nxyz
      integer nyz1,nxy2,nxyz1,nxyz2,n1,n2

      real ta1,ta2,ta3
      common /scrns/ ta1 (lx1*ly1*lz1,lelv)
     $ ,             ta2 (lx1*ly1*lz1,lelv)
     $ ,             ta3 (lx1*ly1*lz1,lelv)

      real           dx   (lx2*ly2*lz2,lelv)
      real           dy   (lx2*ly2*lz2,lelv)
      real           x    (lx1*ly1*lz1,lelv)
      real           rm2  (lx2*ly2*lz2,lelv)
      real           sm2  (lx2*ly2*lz2,lelv)
      real           tm2  (lx2*ly2*lz2,lelv)

      integer e
      real iegr

      nxyz  = lx1*ly1*lz1
      ntot1 = nxyz*nelv
      ntot2 = lx2*ly2*lz2*nelv

      if (istep.eq.0) then

        call outpost(v1mask,v2mask,v3mask,pr,v3mask,'msk')

!         nfaces = 2*ndim
!         do ie=1,nelv
!         do iface=1,nfaces
!            cb  = cbc(iface,ie,ifield)
! !           bc1 = bc(1,iface,ie,ifield)
! !           bc2 = bc(2,iface,ie,ifield)
! !           bc3 = bc(3,iface,ie,ifield)
! 
!            if (cb.ne.'E  ') then
!              write(6,*) ie,iface,cb
!            endif
! 
!         enddo
!         enddo

!        call opcopy(tmp1,tmp2,tmp3,vx,vy,vz)
!        call opdsop(tmp1,tmp2,tmp3,'MXA')    

!        call outpost(vx,vy,vz,pr,t,'  ')
        ifield = 1
!        call rmask(vx,vy,vz,nelv)
!        call outpost(vx,vy,vz,pr,t,'  ')

        if (nio.eq.0) write(6,*) 'NFIELD', nfield
        if (nio.eq.0) write(6,*) 'IFADVC', (ifadvc(i),i=1,nfield)

        ifstrs = .false.
        call bcmask
        ifstrs = .true.

!       Forced to zero. 
        nfaces = 2*ldim
        nedge = 12
        ncrnr = 2**ldim
        ntotf = nface*nelv
        ntots = nedge*nelv
        ntotc = ncrnr*nelv
       
        iflmsf(1) = .false.
        iflmse(1) = .false.
        iflmsc(1) = .false.
        call lfalse (ifmsfc(1,1,1),ntotf)
        call lfalse (ifmseg(1,1,1),ntots)
        call lfalse (ifmscr(1,1,1),ntotc)
      
        call copy(v3mask,v2mask,ntot1)      
        call outpost(v1mask,v2mask,v3mask,pr,v3mask,'msk')

        call col2(vx,v1mask,ntot1)
        call col2(vy,v2mask,ntot1)
        call col2(vz,v3mask,ntot1)
       
        call outpost(vx,vy,vz,pr,vz,'   ')

        call fs_gen_damping


        call rzero(ta1,ntot1)
        call rzero(ta2,ntot1)
        call rzero(ta3,ntot1)
        do i=1,nelv
          iegr = lglel(i)+0.0
          write(6,*) 'ieg', iegr
          call cfill(ta1(1,i),iegr,nxyz)
        enddo
        call outpost(ta1,ta2,ta3,pr,ta3,'eln') 

!        call exitt
      endif  

!      if (istep.eq.1) then

!      endif  


!      if (istep.eq.1) then
!        call outpost(tmp1,tmp2,tmp3,pr,tmp3,'sol')
!        call outpost(tmp5,tmp6,tmp7,pr,tmp7,'sol')
!        call exitt
!      endif  



      return
      end subroutine test_random
!---------------------------------------------------------------------- 

      subroutine writew3m2

      implicit none

      include 'SIZE'
      include 'WZ'

      integer i,j

       do j=1,ly2
         write(6,'(A4,1x,6(F7.5,1x))') 'w3m2', (w3m2(i,j,1), i=1,lx2)
       enddo  

      return
      end subroutine writew3m2
!---------------------------------------------------------------------- 


c automatically added by makenek
      subroutine usrdat0() 

      return
      end

c automatically added by makenek
      subroutine usrsetvert(glo_num,nel,nx,ny,nz) ! to modify glo_num
      integer*8 glo_num(1)

      return
      end

c automatically added by makenek
      subroutine userqtl

      call userqtl_scig

      return
      end
