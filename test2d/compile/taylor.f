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
            
        call gen_mapping_mvb
       
      endif  

      call frame_monitor

      ifto = .true.

      ntot1 = lx1*ly1*lz1*nelv
      ntot2 = lx2*ly2*lz2*lelv


      if (istep.eq.0) then
!        call outpost(tmp1,tmp2,tmp3,pr,tmp3,'ini')
!        call copy(vz,t,ntot1)
       
        call initp_f3d
        ifheat = .false.

      endif

      if (istep.gt.0) then
        call copy(t,vz,ntot1)
      endif  

      call test_constrain()
      call exitt

      call test_random

      if (fs_iffs) call fs_mvmesh()
!      call rzero3(wx,wy,wz,ntot1)

      if (istep.eq.nsteps.or.lastep.eq.1) then
        call frame_end
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine userbc (ix,iy,iz,iside,ieg)

      implicit none

      include 'SIZE'
      include 'NEKUSE'
      include 'GEOM'

      integer ix,iy,iz,iside,ieg
      real pi

      integer jp
      common /ppointr/ jp

      real rmid

      real glmin,glmax
      integer n

      real rad1,rad2,omega1,omega2
      real a1,a2
      common /cylindrical/ rad1,rad2,omega1,omega2,a1,a2


      if (jp.eq.0) then
        ux   = 0.0
        uy   = 0.0
        uz   = 0.0

        rmid = (rad1+rad2)/2
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

      integer isave
      save isave
      data isave /0/

      real glmin,glmax
      integer n

      real rad1,rad2,omega1,omega2
      real a1,a2
      common /cylindrical/ rad1,rad2,omega1,omega2,a1,a2

      real rmid

      ifcouette         = .false.
      ifpoiseuille      = .false.
      iftaylor          = .false.
      iftestmvb         = .true.

      pi = 4.0*atan(1.0)


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
          rmid = (rad1+rad2)/2.0
          ux   = -uparam(1)*exp(-((y-rmid)/0.25)**2)
          uy   = -0.1*uparam(1)*exp(-((y-rmid)/0.25)**2)
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

      implicit none        

      include 'SIZE'
      include 'SOLN'    ! tmult
      include 'INPUT'
      include 'GEOM'

      integer n
      real rad1,rad2,omega1,omega2
      real a1,a2
      common /cylindrical/ rad1,rad2,omega1,omega2,a1,a2

      real glmin,glmax

!      include 'TOTAL'

!      param(54) = -1  ! use >0 for const flowrate or <0 bulk vel
                      ! flow direction is given by (1=x, 2=y, 3=z) 
!      param(55) = 1.0 ! flowrate/bulk-velocity 

!      call gen_mapping_mvb

      n = lx1*ly1*lz1*nelv
      rad1 = glmin(ym1,n)
      rad2 = glmax(ym1,n)
      omega1 = 1.0/rad1
      omega2 = 0.0
      a1 = (omega2*(rad2**2) - omega1*rad1*rad1)/(rad2**2 - rad1**2)
      a2 = (omega1 - omega2)*(rad1**2)*(rad2**2)/(rad2**2 - rad1**2)

      if (nio.eq.0) write(6,*) 'Cylindrical Params:', rad1,rad2,a1,a2

      return
      end
c-----------------------------------------------------------------------
!=======================================================================

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
      include 'MVGEOM'

      include 'FS_ALE'
      include 'SUBGEOM'

      integer ntot1,ntot2
      integer e,i,j

      integer nxyz
      integer nxyz9,ntot9


      real ta1,ta2,ta3
      common /scrns/ ta1 (lx1*ly1*lz1,lelv)
     $ ,             ta2 (lx1*ly1*lz1,lelv)
     $ ,             ta3 (lx1*ly1*lz1,lelv)

      real x,y
      integer idir
      integer igeom

      real om           ! angular frequency


      nxyz  = lx1*ly1*lz1
      ntot1 = nxyz*nelv
      ntot2 = lx2*ly2*lz2*nelv

      nxyz9 = lx9*ly9*lz9
      ntot9 = nxyz9*nelv

      if (istep.eq.0) then

        om = 1.0*pi

        do i=1,ntot1
          x = xm1(i,1,1,1)
          y = ym1(i,1,1,1)
          ym1(i,1,1,1) = ym1(i,1,1,1) + 0.1*cos(om*x) + 0.1*sin(om*y)
          xm1(i,1,1,1) = xm1(i,1,1,1) + 0.05*sin(om*y) + 0.2*cos(om*x) 
        enddo
        
        igeom = 2
        istep = 1
        call gengeom(igeom)

        do i=1,ntot1
          x = xm1(i,1,1,1)
          vx(i,1,1,1) = sxm1(i,1,1,1)
        enddo

        call copy(vx,jacm1,ntot1)

        call subgeom_setup_dssum()

        call genwz_subp()     ! Reference matrices
        call genxyz_subp(xm9,ym9,zm9,lx9,ly9,lz9)

        do i=1,ntot9
          x = xm9(i,1,1,1)
          y = ym9(i,1,1,1)
          ym9(i,1,1,1) = ym9(i,1,1,1) + 0.1*cos(om*x) + 0.1*sin(om*y)
          xm9(i,1,1,1) = xm9(i,1,1,1) + 0.05*sin(om*y) + 0.2*cos(om*x)
        enddo
        call geom9_subp()

!        call fgslib_gs_op(subgeom_gs_handle,xm9,1,1,0)  ! 1 ==> +

        do i=1,ntot9
          x = sxm9(i,1,1,1)
          ta1(i,1) = x
        enddo
        call copy(ta1,jacm9,ntot9)  

        idir = 1        ! 9 -> 1
        call map91(vy,ta1,idir)

        call outpost(vx,vy,vz,pr,t,'tst')

        call exitt
      endif  


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

      subroutine test_constrain()

      implicit none

      include 'SIZE'
      include 'FS_ALE'
      include 'GEOM'
      include 'SOLN'

      include 'MVGEOM'

      integer lt
      parameter (lt=lx1*ly1*lz1*lelt)

      integer i,j,k,e,n,m
      real x,y,z

      real tb1,tb2
      common /scrch/ tb1     (lt)
     $ ,             tb2     (lt)

      real ta1,ta2,ta3,ta4
      common /scrmg/ ta1    (lt)
     $ ,             ta2    (lt)
     $ ,             ta3    (lt)
     $ ,             ta4    (lt)

!     Temporary arrays
      real wk1,wk2,wk3,wk4
      common /scrns/ wk1(lt*2),
     $               wk2(lt*2),
     $               wk3(lt*2),
     $               wk4(lt) 
    
      integer maxiter 
      real vlsc2
      real temp1,temp2

      n = nx1*ny1*nz1*nelv

      do i=1,n
        x           = xm1(i,1,1,1)
        y           = ym1(i,1,1,1)
        wx(i,1,1,1) = cos(y+1)
        wy(i,1,1,1) = sin(y+1)
      enddo


      call fs_global_basis()

      call fs_gllo_xyz

      call fs_gllo_flds(wx,wy,wz)
      call fs_intp_setup
      call fs_get_localpts      

      m = lxfs*lyfs

      call rzero(tb1,m)
      call rzero(tb2,m)

      k = 0
      do j=1,lyfs
      do i=1,lxfs
        k = k+1
        tb1(k)=fld_fs(i,j,1)
      enddo
      enddo
!      call copy(tb1,fld_fs,m)
      call col2(tb1,bm_fs,m)

      k = 0
      do j=1,lyfs
      do i=1,lxfs
        k = k + 1
        write(6,'(A5,2x,6(E16.8E2,2x))') 'glpts',xg_fs(i,j,1),
     $      yg_fs(i,j,1),fld_fs(i,j,1),bm_fs(i,j),
     $      tb1(k),tb2(k)
      enddo
      enddo

      maxiter = 1000
      call fs_gmres(tb1,tb2,maxiter)

      call tensory_op(ta1,tb1,lxfs,lyfs,1,dyt_fs,lyfs)
     
      k = 0
      do j=1,lyfs
      do i=1,lxfs
        k = k + 1
        write(6,'(A6,2x,6(E16.8E2,2x))') 'soln1',xg_fs(i,j,1),
     $      yg_fs(i,j,1),fld_fs(i,j,1),tb1(k),tb2(k),ta1(k)
      enddo
      enddo

      call fs_lagrangian_deriv(ta1,ta2,tb1,tb2,wk1)
      call subcol3(ta1,fld_fs,bm_fs,m)
      temp1 = vlsc2(ta1,ta1,m)
      temp2 = vlsc2(ta2,ta2,m)

      write(6,*) '||Lx-b|| = ', temp1,temp2


      return
      end subroutine
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
