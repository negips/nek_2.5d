!---------------------------------------------------------------------- 
!     Author:      Prabal Negi
!     Description: Routines to generate sub-parametric mesh
!                  using algorithms from the Book
!                  Implementing Spectral Methods for 
!                  Partial Differential Equations: 
!                  Algorithms for Scientists and Engineers (2009)       
!                  - David A. Kopriva
!
!     Note: I have kept things as independent of speclib as I can.
!           If this goes into production runs, I will make it
!           completely independent.      
!
!     Routines:
!     1) genwz_subp():   Generate Reference element matrices      
!      
!---------------------------------------------------------------------- 
      subroutine genwz_subp
C-----------------------------------------------------------------
C
C     GENERATE
C
C            - DERIVATIVE OPERATORS
C            - INTERPOLATION OPERATORS
C            - WEIGHTS
C            - COLLOCATION POINTS
C
C     ASSOCIATED WITH THE
C
C            - GAUSS-LOBATTO LEGENDRE MESH (SUFFIX M1/M2/M3)
C            - GAUSS LEGENDRE         MESH (SUFFIX M2)
C            - GAUSS-LOBATTO JACOBI   MESH (SUFFIX M1/M2/M3)
C
C-----------------------------------------------------------------

      implicit none

      INCLUDE 'SIZE'
      INCLUDE 'WZ'
      INCLUDE 'DXYZ'
      INCLUDE 'IXYZ'
      INCLUDE 'INPUT'
      INCLUDE 'SUBGEOM'

      integer ix,iy,iz

!     Gauss-Lobatto Legendre mesh (suffix M9)
!     Generate collocation points and weights
      call zwgll (zgm9(1,1),wxm9,lx9)
      call zwgll (zgm9(1,2),wym9,ly9)

!     Generating Baryweights for Mesh 1          
      call BarycentricWeights(bwxm1,zgm1(1,1),lx1)
      call BarycentricWeights(bwym1,zgm1(1,2),ly1)

!     Generating Baryweights for Mesh 2 
      call BarycentricWeights(bwxm2,zgm2(1,1),lx2)
      call BarycentricWeights(bwym2,zgm2(1,2),ly2)

!     Generating Baryweights for Mesh 9 
      call BarycentricWeights(bwxm9,zgm9(1,1),lx9)
      call BarycentricWeights(bwym9,zgm9(1,2),ly9)

!     Generate interpolation operators
!     Mesh 9 to Mesh 1        
      call PolynomialInterpolationMatrix(ixm91,zgm1(1,1),lx1,
     $       zgm9(1,1),bwxm9,lx9)
      call transpose(ixtm91,lx9,ixm91,lx1)
      call PolynomialInterpolationMatrix(iym91,zgm1(1,2),ly1,
     $       zgm9(1,2),bwym9,ly9)
      call transpose(iytm91,ly9,iym91,ly1)

!     Mesh 1 to Mesh 9
      call PolynomialInterpolationMatrix(ixm19,zgm9(1,1),lx9,
     $       zgm1(1,1),bwxm1,lx1)
      call transpose(ixtm19,lx1,ixm19,lx9)
      call PolynomialInterpolationMatrix(iym19,zgm9(1,2),ly9,
     $       zgm1(1,2),bwym1,ly1)
      call transpose(iytm19,ly1,iym19,ly9)

!     Generate interpolation operators
!--------------------------------------------------       
!     Mesh 9 to Mesh 2 
      call PolynomialInterpolationMatrix(ixm92,zgm2(1,2),lx2,
     $       zgm9(1,1),bwxm9,lx9)
      call transpose(ixtm92,lx9,ixm92,lx2)
      call PolynomialInterpolationMatrix(iym92,zgm2(1,2),ly2,
     $       zgm9(1,2),bwym9,ly9)
      call transpose(iytm92,ly9,iym92,ly2)
      izm91(lz2,lz9)  = 1.0
      iztm91(lz9,lz2) = 1.0

!     Mesh 2 to Mesh 9
      call PolynomialInterpolationMatrix(ixm29,zgm9(1,1),lx9,
     $       zgm2(1,1),bwxm2,lx2)
      call transpose(ixtm29,lx2,ixm29,lx9)
      call PolynomialInterpolationMatrix(iym29,zgm9(1,2),ly9,
     $       zgm2(1,2),bwym2,ly2)
      call transpose(iytm29,ly2,iym29,ly9)

!     Derivative Matrices
!--------------------------------------------------         
!     Derivative matrices at the Nodal points
      call PolynomialDerivativeMatrix(dxm9,zgm9(1,1),bwxm9,lx9)
      call transpose(dxtm9,lx9,dxm9,lx9)
      call PolynomialDerivativeMatrix(dym9,zgm9(1,2),bwym9,ly9)
      call transpose(dytm9,ly9,dym9,ly9)

!     Derivative operators for the staggered mesh
!     Mesh 9 to Mesh 1        
      call LagrangeDerivativeMatrix(dxm91,zgm1(1,1),lx1,
     $                              zgm9(1,1),bwxm9,lx9)
      call transpose(dxtm91,lx9,dxm91,lx1)
      call LagrangeDerivativeMatrix(dym91,zgm1(1,2),ly1,
     $                              zgm9(1,2),bwym9,ly9)
      call transpose(dytm91,ly9,dym91,ly1)

!     Mesh 1 to Mesh 9
      call LagrangeDerivativeMatrix(dxm19,zgm9(1,1),lx9,
     $                              zgm1(1,1),bwxm1,lx1)
      call transpose(dxtm19,lx1,dxm19,lx9)
      call LagrangeDerivativeMatrix(dym19,zgm9(1,2),ly9,
     $                              zgm1(1,2),bwym1,ly1)
      call transpose(dytm19,ly1,dym19,ly9)

!     Mesh 9 to Mesh 2 
      call LagrangeDerivativeMatrix(dxm92,zgm2(1,1),lx2,
     $                              zgm9(1,1),bwxm9,lx9)
      call transpose(dxtm92,lx9,dxm92,lx2)
      call LagrangeDerivativeMatrix(dym92,zgm2(1,2),ly2,
     $                              zgm9(1,2),bwym9,ly9)
      call transpose(dytm92,ly9,dym92,ly2)

!     Mesh 2 to Mesh 9
      call LagrangeDerivativeMatrix(dxm29,zgm9(1,1),lx9,
     $                              zgm2(1,1),bwxm2,lx2)
      call transpose(dxtm29,lx2,dxm29,lx9)
      call LagrangeDerivativeMatrix(dym29,zgm9(1,2),ly9,
     $                              zgm2(1,2),bwym2,ly2)
      call transpose(dytm29,ly2,dym29,ly9)

      IF (ldim.EQ.2) THEN

!       Two-dimensional case 

!       Gauss-Lobatto Legendre mesh (suffix M9)
        zgm9(lz9,3) = 0.
        wzm9(lz9)   = 1.
        do iy=1,ly9
        do ix=1,lx9
          w3m9(ix,iy,1)=wxm9(ix)*wym9(iy)
        enddo
        enddo  

!       Generating Baryweights for Mesh 1          
        bwzm1(lz1) = 1.0

!       Generating Baryweights for Mesh 2 
        bwzm2(lz2) = 1.0

!       Generating Baryweights for Mesh 9 
        bwzm9(lz9) = 1.0

!       Generate interpolation operators
!       Mesh 9 to Mesh 1        
        izm91(lz1,lz9)  = 1.0
        iztm91(lz9,lz1) = 1.0

!       Mesh 1 to Mesh 9
        izm19(lz9,lz1)  = 1.0
        iztm19(lz1,lz9) = 1.0

!       Generate interpolation operators
!       Mesh 9 to Mesh 2 
        izm91(lz2,lz9)  = 1.0
        iztm91(lz9,lz2) = 1.0

!       Mesh 2 to Mesh 9
        izm19(lz9,lz2)  = 1.0
        iztm19(lz2,lz9) = 1.0

!       Derivative Matrices
!--------------------------------------------------         

!       Compute derivative matrices at the Nodal points
        call rzero (dzm9 ,lz9*lz9)
        call rzero (dztm9,lz9*lz9)

!       Compute derivative operators for the staggered mesh
!       Mesh 9 to Mesh 1        
        call rzero (dzm91 ,lz1*lz9)
        call rzero (dztm91,lz9*lz1)

!       Mesh 1 to Mesh 9
        call rzero (dzm19 ,lz9*lz1)
        call rzero (dztm19,lz1*lz9)

!       Mesh 9 to Mesh 2 
        call rzero (dzm92 ,lz2*lz9)
        call rzero (dztm92,lz9*lz2)

!       Mesh 2 to Mesh 9
        call rzero (dzm29 ,lz9*lz2)
        call rzero (dztm29,lz2*lz9)

      ELSE
!       Three-dimensional case

!       Gauss-Lobatto Legendre mesh (suffix M9)
!       Generate collocation points and weights
        call zwgll (zgm9(1,3),wzm9,lz9)
        do iz=1,lz9
        do iy=1,ly9
        do ix=1,lx9
          w3m9(ix,iy,1)=wxm9(ix)*wym9(iy)*wzm9(iz)
        enddo
        enddo
        enddo

!       Generating Baryweights for Mesh 1          
        call BarycentricWeights(bwzm1,zgm1(1,3),lz1)

!       Generating Baryweights for Mesh 2 
        call BarycentricWeights(bwzm2,zgm2(1,3),lz2)

!       Generating Baryweights for Mesh 9 
        call BarycentricWeights(bwzm9,zgm9(1,3),lz9)

!       Generate interpolation operators
!       Mesh 9 to Mesh 1        
        call PolynomialInterpolationMatrix(izm91,zgm1(1,3),lz1,
     $         zgm9(1,3),bwzm9,lz9)
        call transpose(iztm91,lz9,izm91,lz1)

!       Mesh 1 to Mesh 9
        call PolynomialInterpolationMatrix(izm19,zgm9(1,3),lz9,
     $         zgm1(1,3),bwzm1,lz1)
        call transpose(iztm19,lz1,izm19,lz9)

!       Generate interpolation operators
!       Mesh 9 to Mesh 2 
        call PolynomialInterpolationMatrix(izm92,zgm2(1,3),lz2,
     $         zgm9(1,3),bwzm9,lz9)
        call transpose(iztm92,lz9,izm92,lz2)

!       Mesh 2 to Mesh 9
        call PolynomialInterpolationMatrix(izm29,zgm9(1,3),lz9,
     $         zgm2(1,3),bwzm2,lz2)
        call transpose(iztm29,lz2,izm29,lz9)

!       Derivative Matrices
!--------------------------------------------------         

!       Compute derivative matrices at the Nodal points
        call PolynomialDerivativeMatrix(dzm9,zgm9(1,3),bwzm9,lz9)
        call transpose(dztm9,lz9,dzm9,lz9)

!       Compute derivative operators for the staggered mesh
!       Mesh 9 to Mesh 1        
        call LagrangeDerivativeMatrix(dzm91,zgm1(1,3),lz1,
     $                                zgm9(1,3),bwzm9,lz9)
        call transpose(dztm91,lz9,dzm91,lz1)

!       Mesh 1 to Mesh 9
        call LagrangeDerivativeMatrix(dzm19,zgm9(1,3),lz9,
     $                                zgm1(1,3),bwzm1,lz1)
        call transpose(dztm19,lz1,dzm19,lz9)

!       Mesh 9 to Mesh 2 
        call LagrangeDerivativeMatrix(dzm92,zgm2(1,3),lz2,
     $                                zgm9(1,3),bwzm9,lz9)
        call transpose(dztm92,lz9,dzm92,lz2)

!       Mesh 2 to Mesh 9
        call LagrangeDerivativeMatrix(dzm29,zgm9(1,3),lz9,
     $                                zgm2(1,3),bwzm2,lz2)
        call transpose(dztm29,lz2,dzm29,lz9)

      endif

      return
      end subroutine
!---------------------------------------------------------------------- 
      subroutine subgeom_setup_dssum

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'

      include 'SUBGEOM'

      integer lxyz
      parameter (lxyz=(lx1+2)*(ly1+2)*(lz1+2))
      common /c_is1/ glo_num(lxyz*lelv)
      common /ivrtx/ vertex ((2**ldim)*lelt)

      integer*8 glo_num
      integer vertex
      integer nx,ny,nz
      
!     set up direct stiffness summation for sub parametric geometry
      call get_vert

      nx=lx9
      ny=ly9
      nz=lz9
      call setupds(subgeom_gs_handle,nx,ny,nz
     $             ,nelv,nelgv,vertex,glo_num)

      return
      end subroutine
c----------------------------------------------------------------------
     
      subroutine gengeom_subp (igeom)

!     Generate geometry data

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
      include 'GEOM'
      include 'WZ'
      include 'SUBGEOM'

      integer igeom

!      if (nio.eq.0.and.istep.le.1) write(6,*) 'generate geometry data'

      IF (IGEOM.EQ.1) THEN
         RETURN
      ELSEIF (IGEOM.EQ.2) THEN
         CALL LAGMASS
         IF (ISTEP.EQ.0) call genxyz_subp(xm9,ym9,zm9,lx9,ly9,lz9)
         IF (ISTEP.GE.1) CALL UPDCOOR
!        Calculate Geometric Matrices
!        dr/dx etc. on M9, JacM9 
!        dx/dr etc. on M1 (Saved in scratch. Used in subsequent routines)
!        dr/dx etc. on M1, JacM1, BM1, unx, t1x, t2x etc, GiMi
         call geom9_subp 
         CALL GEOM2           ! Geometric factors on Mesh 2
         CALL UPDMSYS (1)     ! Update Masks
         CALL VOLUME          ! Calculate volume
         CALL SETINVM         ! Inverse Mass Matrices 
         CALL SETDEF          ! Logicals for deformed elements 
         CALL SFASTAX         ! Set Matrices for undeformed elements
      ENDIF

      if (nio.eq.0.and.istep.le.1) then
        write(6,*) 'done :: generate geometry data' 
        write(6,*) ' '
      endif

      return
      end
c-----------------------------------------------------------------------

      subroutine genxyz_subp(xml,yml,zml,nxl,nyl,nzl)

      implicit none  

      include 'SIZE'
      include 'WZ'
      include 'GEOM'
      include 'TOPOL'
      include 'INPUT'
      include 'PARALLEL'

      integer nxl,nyl,nzl
      real xml(nxl,nyl,nzl,1),yml(nxl,nyl,nzl,1),zml(nxl,nyl,nzl,1)

!     Note : CTMP1 is used in this format in several subsequent routines
      real h,xcrved,ycrved,zcrved,zgml,work
      common /ctmp1/ h(lx1,3,2),xcrved(lx1),ycrved(ly1),zcrved(lz1)
     $             , zgml(lx1,3),work(3,lx1,lz1)

      integer ldw
      parameter (ldw=2*lx1*ly1*lz1)
      real w
      common /ctmp0/ w(ldw)

      integer ie,isid,iface,nfaces

      character*1 ccv

c     Initialize geometry arrays with bi- triquadratic deformations
      call linquad_subp(xml,yml,zml,nxl,nyl,nzl)

      do ie=1,nelt

         call setzgml_subp(zgml,ie,nxl,nyl,nzl,ifaxis)
         call sethmat (h,zgml,nxl,nyl,nzl)

c        Deform surfaces - general 3D deformations
c                        - extruded geometry deformations
         nfaces = 2*ldim
         do iface=1,nfaces
           ccv = ccurve(iface,ie)
           if (ccv.eq.'s') 
     $        call sphsrf(xml,yml,zml,iface,ie,nxl,nyl,nzl,work) 
           if (ccv.eq.'e') 
     $        call gensrf(xml,yml,zml,iface,ie,nxl,nyl,nzl,zgml) 
         enddo

         do isid=1,8
           ccv = ccurve(isid,ie)
           if (ccv.eq.'C') call arcsrf(xml,yml,zml,nxl,nyl,nzl,ie,isid)
         enddo

      enddo

      return
      end
c-----------------------------------------------------------------------

      subroutine setzgml_subp (zgml,e,nxl,nyl,nzl,ifaxl)

      include 'SIZE'
      include 'WZ'
      include 'GEOM'
      include 'SUBGEOM'

      real zgml(lx1,3)
      integer e
      logical ifaxl

      call rzero (zgml,3*lx1)

      if (nxl.eq.3 .and. .not. ifaxl) then
         do k=1,3
            zgml(1,k) = -1
            zgml(2,k) =  0
            zgml(3,k) =  1
         enddo
      elseif (ifgmsh3.and.nxl.eq.lx3) then
         call copy(zgml(1,1),zgm3(1,1),lx3)
         call copy(zgml(1,2),zgm3(1,2),ly3)
         call copy(zgml(1,3),zgm3(1,3),lz3)
         if (ifaxl .and. ifrzer(e)) call copy(zgml(1,2),zam3,ly3)
      elseif (nxl.eq.lx1) then
         call copy(zgml(1,1),zgm1(1,1),lx1)
         call copy(zgml(1,2),zgm1(1,2),ly1)
         call copy(zgml(1,3),zgm1(1,3),lz1)
         if (ifaxl .and. ifrzer(e)) call copy(zgml(1,2),zam1,ly1)
      else
         call copy(zgml(1,1),zgm9(1,1),lx9)
         call copy(zgml(1,2),zgm9(1,2),ly9)
         call copy(zgml(1,3),zgm9(1,3),lz9)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine linquad_subp(xl,yl,zl,nxl,nyl,nzl)

      implicit none  

      include 'SIZE'
      include 'WZ'
      include 'GEOM'
      include 'TOPOL'
      include 'INPUT'
      include 'PARALLEL'

      integer nxl,nyl,nzl
      real xl(nxl*nyl*nzl,1),yl(nxl*nyl*nzl,1),zl(nxl*nyl*nzl,1)

      integer e,k
      logical ifmid

      integer nedge

      nedge = 4 + 8*(ldim-2)

      do e=1,nelt ! Loop over all elements

         ifmid = .false.
         do k=1,nedge
            if (ccurve(k,e).eq.'m') ifmid = .true.
         enddo

         if (lx1.eq.2) ifmid = .false.
         if (ifmid) then
            call xyzquad_subp(xl(1,e),yl(1,e),zl(1,e),nxl,nyl,nzl,e)
         else
            call xyzlin_subp (xl(1,e),yl(1,e),zl(1,e),
     $                        nxl,nyl,nzl,e,ifaxis)
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------

      subroutine xyzlin_subp(xl,yl,zl,nxl,nyl,nzl,e,ifaxl)
c     Generate bi- or trilinear mesh

      implicit none

      include 'SIZE'
      include 'INPUT'

      integer nxl,nyl,nzl
      real xl(nxl,nyl,nzl),yl(nxl,nyl,nzl),zl(nxl,nyl,nzl)
      integer e,i,k,ix,ldim2
      logical ifaxl ! local ifaxis specification

c   Preprocessor Corner notation:      Symmetric Corner notation:
c
c           4+-----+3    ^ s                    3+-----+4    ^ s
c           /     /|     |                      /     /|     |
c          /     / |     |                     /     / |     |
c        8+-----+7 +2    +----> r            7+-----+8 +2    +----> r
c         |     | /     /                     |     | /     /
c         |     |/     /                      |     |/     /
c        5+-----+6    t                      5+-----+6    t

      integer indx(8)
      save    indx
      data    indx / 1,2,4,3,5,6,8,7 /

      integer ldw
      parameter (ldw=4*lx1*ly1*lz1)
      real xcb,ycb,zcb,w
      common /ctmp0/ xcb(2,2,2),ycb(2,2,2),zcb(2,2,2),w(ldw)

      real zgml,jx,jy,jz,jxt,jyt,jzt,zlin
      common /cxyzl/ zgml(lx1,3),jx (lx1*2),jy (lx1*2),jz (lx1*2)
     $                          ,jxt(lx1*2),jyt(lx1*2),jzt(lx1*2)
     $                          ,zlin(2)

      call setzgml_subp (zgml,e,nxl,nyl,nzl,ifaxl)

      zlin(1) = -1
      zlin(2) =  1

      k = 1
      do i=1,nxl
         call fd_weights_full(zgml(i,1),zlin,1,0,jxt(k))
         call fd_weights_full(zgml(i,2),zlin,1,0,jyt(k))
         call fd_weights_full(zgml(i,3),zlin,1,0,jzt(k))
         k=k+2
      enddo
      call transpose(jx,nxl,jxt,2)

      ldim2 = 2**ldim
      do ix=1,ldim2          ! Convert prex notation to lexicographical
         i=indx(ix)
         xcb(ix,1,1)=xc(i,e)
         ycb(ix,1,1)=yc(i,e)
         zcb(ix,1,1)=zc(i,e)
      enddo

c     Map R-S-T space into physical X-Y-Z space.

      ! NOTE:  Assumes nxl=nyl=nzl !

      call tensr3(xl,nxl,xcb,2,jx,jyt,jzt,w)
      call tensr3(yl,nxl,ycb,2,jx,jyt,jzt,w)
      call tensr3(zl,nxl,zcb,2,jx,jyt,jzt,w)

      return
      end
c-----------------------------------------------------------------------

      subroutine xyzquad_subp(xl,yl,zl,nxl,nyl,nzl,e)
c     Generate bi- or trilinear mesh

      implicit none

      include 'SIZE'
      include 'INPUT'

      integer nxl,nyl,nzl
      real xl(nxl,nyl,nzl),yl(nxl,nyl,nzl),zl(nxl,nyl,nzl)
      real xq(27),yq(27),zq(27)
      integer e,i,j,k,ix,ldim2,nedge

      integer ldw
      parameter (ldw=4*lx1*ly1*lz1)
      real w,zg
      common /ctmp0/ w(ldw,2),zg(3)

c     Note : CTMP1 is used in this format in several subsequent routines

      integer eindx(12)  ! index of 12 edges into 3x3x3 tensor
      save    eindx      ! Follows preprocessor notation..
      data    eindx /  2 ,  6 ,  8 ,  4
     $              , 20 , 24 , 26 , 22
     $              , 10 , 12 , 18 , 16  /  ! preproc. vtx notation

      real zgml,jx,jy,jz,jxt,jyt,jzt,zquad
      common /cxyzl/ zgml(lx1,3),jx (lx1*3),jy (lx1*3),jz (lx1*3)
     $                          ,jxt(lx1*3),jyt(lx1*3),jzt(lx1*3)
     $                          ,zquad(3)

      call xyzlin(xq,yq,zq,3,3,3,e,.false.) ! map bilin to 3x3x3

      nedge = 4 + 8*(ldim-2)

      do k=1,nedge
         if (ccurve(k,e).eq.'m') then
            j = eindx(k)
            xq(j) = curve(1,k,e)
            yq(j) = curve(2,k,e)
            zq(j) = curve(3,k,e)
         endif
      enddo

      zg(1) = -1
      zg(2) =  0 
      zg(3) =  1

      if (if3d) then
         call gh_face_extend(xq,zg,3,2,w(1,1),w(1,2)) ! 2 --> edge extend
         call gh_face_extend(yq,zg,3,2,w(1,1),w(1,2))
         call gh_face_extend(zq,zg,3,2,w(1,1),w(1,2))
      else
         call gh_face_extend_2d(xq,zg,3,2,w(1,1),w(1,2)) ! 2 --> edge extend
         call gh_face_extend_2d(yq,zg,3,2,w(1,1),w(1,2))
      endif
      call clean_xyzq(xq,yq,zq,if3d) ! verify that midside node is in "middle"


c     Map R-S-T space into physical X-Y-Z space.
      ! NOTE:  Assumes nxl=nyl=nzl !

      zquad(1) = -1
      zquad(2) =  0
      zquad(3) =  1

      call setzgml (zgml,e,nxl,nyl,nzl,ifaxis)  ! Here we address axisymm.

      k = 1
      do i=1,nxl
         call fd_weights_full(zgml(i,1),zquad,2,0,jxt(k))
         call fd_weights_full(zgml(i,2),zquad,2,0,jyt(k))
         call fd_weights_full(zgml(i,3),zquad,2,0,jzt(k))
         k=k+3
      enddo
      call transpose(jx,nxl,jxt,3)

      call tensr3(xl,nxl,xq,3,jx,jyt,jzt,w)
      call tensr3(yl,nxl,yq,3,jx,jyt,jzt,w)
      call tensr3(zl,nxl,zq,3,jx,jyt,jzt,w)

      return
      end
c-----------------------------------------------------------------------
      subroutine geom9_subp()

!     Routine to generate all elemental geometric data for mesh 1.

!     Velocity formulation : global-to-local mapping based on mesh 3
!     Stress   formulation : global-to-local mapping based on mesh 1

      implicit none        

      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
      include 'SUBGEOM'

!     This fills up scrns and ctmp1 scratch arrays
      call glmapm9

      call geodat91

      return
      end subroutine
!---------------------------------------------------------------------- 
      subroutine glmapm9

!     Routine to generate mapping data based on mesh 9
!     (Gauss-Legendre Lobatto meshes).
!
!     XRM9,  YRM9,  ZRM9   -   dx/dr, dy/dr, dz/dr
!     XSM9,  YSM9,  ZSM9   -   dx/ds, dy/ds, dz/ds
!     XTM9,  YTM9,  ZTM9   -   dx/dt, dy/dt, dz/dt
!     RXM9,  RYM9,  RZM9   -   dr/dx, dr/dy, dr/dz
!     SXM9,  SYM9,  SZM9   -   ds/dx, ds/dy, ds/dz
!     TXM9,  TYM9,  TZM9   -   dt/dx, dt/dy, dt/dz
!     JACM1                -   Jacobian

      implicit none

      INCLUDE 'SIZE'
      INCLUDE 'GEOM'
      INCLUDE 'INPUT'
      INCLUDE 'SOLN'
      INCLUDE 'SUBGEOM'

C     Note: Subroutines GLMAPM1, GEODAT1, AREA2, SETWGTR and AREA3 
C           share the same array structure in Scratch Common /SCRNS/.

      REAL XRM9(LX9,LY9,LZ9)
     $ ,   YRM9(LX9,LY9,LZ9)
     $ ,   XSM9(LX9,LY9,LZ9)
     $ ,   YSM9(LX9,LY9,LZ9)
     $ ,   XTM9(LX9,LY9,LZ9)
     $ ,   YTM9(LX9,LY9,LZ9)
     $ ,   ZRM9(LX9,LY9,LZ9)
     $ ,   ZSM9(LX9,LY9,LZ9)
     $ ,   ZTM9(LX9,LY9,LZ9)

      REAL XRM1,YRM1,XSM1,YSM1,XTM1,YTM1,ZRM1,ZSM1,ZTM1
      COMMON /SCRNS/ XRM1(LX1,LY1,LZ1,LELT)
     $ ,             YRM1(LX1,LY1,LZ1,LELT)
     $ ,             XSM1(LX1,LY1,LZ1,LELT)
     $ ,             YSM1(LX1,LY1,LZ1,LELT)
     $ ,             XTM1(LX1,LY1,LZ1,LELT)
     $ ,             YTM1(LX1,LY1,LZ1,LELT)
     $ ,             ZRM1(LX1,LY1,LZ1,LELT)
      COMMON /CTMP1/ ZSM1(LX1,LY1,LZ1,LELT)
     $ ,             ZTM1(LX1,LY1,LZ1,LELT)

      integer e

      integer ntot1
      integer nxy9,nyz9,nxyz9,ntot9
      integer kerr,ierr
      integer iglsum

      integer idir

      nxy9  = lx9*ly9
      nyz9  = ly9*lz9
      nxyz9 = lx9*ly9*lz9
      ntot9 = nxyz9*nelt

      ntot1 = lx1*ly1*lz1*nelt

      do e=1,nelt

!       single element dx/dr etc.          
        call xyzrst9_e(xrm9,yrm9,zrm9,xsm9,ysm9,zsm9,xtm9,ytm9,ztm9,
     $             ifaxis,e)

        if (ldim.eq.2) then
          call rzero  (jacm9(1,1,1,e),             nxyz9)
          call addcol3(jacm9(1,1,1,e), xrm9, ysm9, nxyz9)
          call subcol3(jacm9(1,1,1,e), xsm9, yrm9, nxyz9)
          call copy   (rxm9 (1,1,1,e), ysm9,       nxyz9)
          call copy   (rym9 (1,1,1,e), xsm9,       nxyz9)
          call chsign (rym9 (1,1,1,e),             nxyz9)
          call copy   (sxm9 (1,1,1,e), yrm9,       nxyz9)
          call chsign (sxm9 (1,1,1,e),             nxyz9)
          call copy   (sym9 (1,1,1,e), xrm9,       nxyz9)
          call rzero  (rzm9 (1,1,1,e),             nxyz9)
          call rzero  (szm9 (1,1,1,e),             nxyz9)
          call rone   (tzm9 (1,1,1,e),             nxyz9)
        else
          call rzero  (jacm9(1,1,1,e),                         nxyz9)
          call addcol4(jacm9(1,1,1,e), xrm9, ysm9, ztm9,       nxyz9)
          call addcol4(jacm9(1,1,1,e), xtm9, yrm9, zsm9,       nxyz9)
          call addcol4(jacm9(1,1,1,e), xsm9, ytm9, zrm9,       nxyz9)
          call subcol4(jacm9(1,1,1,e), xrm9, ytm9, zsm9,       nxyz9)
          call subcol4(jacm9(1,1,1,e), xsm9, yrm9, ztm9,       nxyz9)
          call subcol4(jacm9(1,1,1,e), xtm9, ysm9, zrm9,       nxyz9)
          call ascol5 (rxm9 (1,1,1,e), ysm9, ztm9, ytm9, zsm9, nxyz9)
          call ascol5 (rym9 (1,1,1,e), xtm9, zsm9, xsm9, ztm9, nxyz9)
          call ascol5 (rzm9 (1,1,1,e), xsm9, ytm9, xtm9, ysm9, nxyz9)
          call ascol5 (sxm9 (1,1,1,e), ytm9, zrm9, yrm9, ztm9, nxyz9)
          call ascol5 (sym9 (1,1,1,e), xrm9, ztm9, xtm9, zrm9, nxyz9)
          call ascol5 (szm9 (1,1,1,e), xtm9, yrm9, xrm9, ytm9, nxyz9)
          call ascol5 (txm9 (1,1,1,e), yrm9, zsm9, ysm9, zrm9, nxyz9)
          call ascol5 (tym9 (1,1,1,e), xsm9, zrm9, xrm9, zsm9, nxyz9)
          call ascol5 (tzm9 (1,1,1,e), xrm9, ysm9, xsm9, yrm9, nxyz9)
        endif

        idir = 1  ! Mesh 9 -> Mesh 1
        call map91_e(xrm1(1,1,1,e), xrm9, idir)    ! xrm1
        call map91_e(xsm1(1,1,1,e), xsm9, idir)    ! xsm1
        call map91_e(xtm1(1,1,1,e), xtm9, idir)    ! xtm1
        call map91_e(yrm1(1,1,1,e), yrm9, idir)    ! yrm1
        call map91_e(ysm1(1,1,1,e), ysm9, idir)    ! ysm1
        call map91_e(ytm1(1,1,1,e), ytm9, idir)    ! ytm1
        call map91_e(zrm1(1,1,1,e), zrm9, idir)    ! zrm1
        call map91_e(zsm1(1,1,1,e), zsm9, idir)    ! zsm1
        call map91_e(ztm1(1,1,1,e), ztm9, idir)    ! ztm1

        call map91_e(rxm1(1,1,1,e), rxm9(1,1,1,e), idir)    ! rxm1
        call map91_e(sxm1(1,1,1,e), sxm9(1,1,1,e), idir)    ! sxm1
        call map91_e(txm1(1,1,1,e), txm9(1,1,1,e), idir)    ! txm1
        call map91_e(rym1(1,1,1,e), rym9(1,1,1,e), idir)    ! rym1
        call map91_e(sym1(1,1,1,e), sym9(1,1,1,e), idir)    ! sym1
        call map91_e(tym1(1,1,1,e), tym9(1,1,1,e), idir)    ! tym1
        call map91_e(rzm1(1,1,1,e), rzm9(1,1,1,e), idir)    ! rzm1
        call map91_e(szm1(1,1,1,e), szm9(1,1,1,e), idir)    ! szm1
        call map91_e(tzm1(1,1,1,e), tzm9(1,1,1,e), idir)    ! tzm1

        call map91_e(jacm1(1,1,1,e), jacm9(1,1,1,e), idir)  ! jacm1
     
      enddo       ! e=1,nelt 

      kerr = 0
      do e=1,nelt
         call chkjac(jacm9(1,1,1,e),nxyz9,e,xm9(1,1,1,e),
     $               ym9(1,1,1,e),zm9(1,1,1,e),ldim,ierr)
         if (ierr.ne.0) kerr = kerr+1
      enddo

      kerr = iglsum(kerr,1)
      if (kerr.gt.0) then
         ifxyo = .true.
         ifvo  = .false.
         ifpo  = .false.
         ifto  = .true.
!        Interpolate Jacobian to velocity mesh            
         call map91(t,jacm9,1) 
         call outpost(vx,vy,vz,pr,t,'xyz')
         if (nid.eq.0) write(6,*) 'jac error'
         if (nid.eq.0) write(6,*) 'jac output in temperature field'
         call exitt
      endif

      call invers2(jacmi9,jacm9,ntot9)
      call invers2(jacmi,jacm1,ntot1)

      return
      end subroutine
!----------------------------------------------------------------------

      subroutine xyzrst9(xrm9,yrm9,zrm9,xsm9,ysm9,zsm9,
     $                     xtm9,ytm9,ztm9,ifaxis)

!     Compute global-to-local derivatives on mesh 1.

      implicit none

      INCLUDE 'SIZE'
      INCLUDE 'GEOM'
      INCLUDE 'SUBGEOM'

      real  xrm9(lx9,ly9,lz9,1),yrm9(lx9,ly9,lz9,1)
     $    , zrm9(lx9,ly9,lz9,1),xsm9(lx9,ly9,lz9,1)
     $    , ysm9(lx9,ly9,lz9,1),zsm9(lx9,ly9,lz9,1)
     $    , xtm9(lx9,ly9,lz9,1),ytm9(lx9,ly9,lz9,1)
     $    , ztm9(lx9,ly9,lz9,1)

      integer e
      logical ifaxis

      do e=1,nelt

        call xyzrst9_e(xrm9(1,1,1,e),yrm9(1,1,1,e),zrm9(1,1,1,e)
     $                ,xsm9(1,1,1,e),ysm9(1,1,1,e),zsm9(1,1,1,e)
     $                ,xtm9(1,1,1,e),ytm9(1,1,1,e),ztm9(1,1,1,e)
     $                ,ifaxis,e)

      enddo

      return
      end subroutine
!---------------------------------------------------------------------- 
      subroutine xyzrst9_e(xrm9,yrm9,zrm9,xsm9,ysm9,zsm9,
     $                     xtm9,ytm9,ztm9,ifaxis,e)

!     Compute global-to-local derivatives on mesh 1
!     for a single element.
!     I assume we are on mesh 9. But it could be made more general
!     if we pass xm9,ym9,zm9, lx,ly,lz and dxm1 etc.       

      implicit none

      INCLUDE 'SIZE'
      INCLUDE 'GEOM'
      INCLUDE 'SUBGEOM'

      real  xrm9(lx9,ly9,lz9), yrm9(lx9,ly9,lz9)
     $    , zrm9(lx9,ly9,lz9), xsm9(lx9,ly9,lz9)
     $    , ysm9(lx9,ly9,lz9), zsm9(lx9,ly9,lz9)
     $    , xtm9(lx9,ly9,lz9), ytm9(lx9,ly9,lz9)
     $    , ztm9(lx9,ly9,lz9)

      logical ifaxis
      integer e         ! element number
      integer nxy9

      nxy9=lx9*ly9

!      if (ifaxis) call setaxdy ( ifrzer(iel) )

!     d()/dr
      call tensorx_op(xrm9,xm9(1,1,1,e),lx9,ly9,lz9,
     $                dxm9,lx9)
      call tensorx_op(yrm9,ym9(1,1,1,e),lx9,ly9,lz9,
     $                dxm9,lx9)
      call tensorx_op(zrm9,zm9(1,1,1,e),lx9,ly9,lz9,
     $                dxm9,lx9)

!     d()/ds
      call tensory_op(xsm9,xm9(1,1,1,e),lx9,ly9,lz9,
     $                dytm9,ly9)
      call tensory_op(ysm9,ym9(1,1,1,e),lx9,ly9,lz9,
     $                dytm9,ly9)
      call tensory_op(zsm9,zm9(1,1,1,e),lx9,ly9,lz9,
     $                dytm9,ly9)

!     d()/dt
      if (ldim.eq.3) then
        call tensorz_op(xtm9,xm9(1,1,1,e),lx9,ly9,lz9,
     $                  dztm9,lz9)
        call tensorz_op(ytm9,ym9(1,1,1,e),lx9,ly9,lz9,
     $                  dztm9,lz9)
        call tensorz_op(ztm9,zm9(1,1,1,e),lx9,ly9,lz9,
     $                  dztm9,lz9)
      else  
         call rzero (xtm9,nxy9)
         call rzero (ytm9,nxy9)
         call rone  (ztm9,nxy9)
      endif

      return
      end subroutine
!---------------------------------------------------------------------- 
      subroutine map91(fld1,fld9,idir)

      implicit none

      include 'SIZE'
      include 'SUBGEOM'

      real fld1(lx1,ly1,lz1,nelt)
      real fld9(lx9,ly9,lz9,nelt)

      integer idir      ! idir: 1 => Forwards; idir: -1 => Backwards;
      integer e

      do e=1,nelt
        call map91_e(fld1(1,1,1,e),fld9(1,1,1,e),idir)
      enddo  

      end subroutine      
!---------------------------------------------------------------------- 
      subroutine map91_e(fld1,fld9,idir)

      implicit none

      include 'SIZE'
      include 'SUBGEOM'

      real fld1(lx1,ly1,lz1)
      real fld9(lx9,ly9,lz9)

      integer idir      ! idir: 1 => Forwards; idir: -1 => Backwards;

!     work array for tensor_op_ifi      
      real tnsr_wk(lx1*ly1*lz1,2)
      common /tnsr_op_work/ tnsr_wk

      logical ifx,ify,ifz

!     If Identity logicals
      ifx = .false.
      ify = .false.
      ifz = .false.
      if (lx1.eq.lx9) ifx = .true.
      if (ly1.eq.ly9) ify = .true.
      if (lz1.eq.lz9) ifz = .true.

      if (idir.eq.1) then      
        call tensor3_op_ifi(fld1,fld9,lx9,ly9,lz9,ixm91
     $      ,iytm91,iztm91,lx1,ly1,lz1,ifx,ify,ifz,tnsr_wk,lx1)
      elseif (idir.eq.-1) then
        call tensor3_op_ifi(fld9,fld1,lx1,ly1,lz1,ixm19
     $      ,iytm19,iztm19,lx9,ly9,lz9,ifx,ify,ifz,tnsr_wk,lx1)
      else
        if (nid.eq.0) write(6,*) 'Unknown idir in map91.'
        call exitt  
      endif  


      return
      end subroutine
!----------------------------------------------------------------------
      subroutine map92(fld2,fld9,idir)

      implicit none

      include 'SIZE'
      include 'SUBGEOM'

      real fld2(lx2,ly2,lz2,nelt)
      real fld9(lx9,ly9,lz9,nelt)

      integer idir      ! idir: 1 => Forwards; idir: -1 => Backwards;
      integer e

      do e=1,nelt
        call map92_e(fld2(1,1,1,e),fld9(1,1,1,e),idir)
      enddo  

      end subroutine      
!---------------------------------------------------------------------- 

      subroutine map92_e(fld2,fld9,idir)

      implicit none

      include 'SIZE'
      include 'SUBGEOM'

      real fld2(lx2,ly2,lz2)
      real fld9(lx9,ly9,lz9)

      integer idir      ! idir: 1 => Forwards; idir: -1 => Backwards;

!     work array for tensor_op_ifi      
      real tnsr_wk(lx1*ly1*lz1,2)
      common /tnsr_op_work/ tnsr_wk

      logical ifx,ify,ifz

!     If Identity logicals
      ifx = .false.
      ify = .false.
      ifz = .false.
!     Always false since Mesh 2 is on Gauss Legendre Mesh.      
!      if (lx2.eq.lx9) ifx = .true.
!      if (ly2.eq.ly9) ify = .true.
!      if (lz2.eq.lz9) ifz = .true.

      if (idir.eq.1) then      
        call tensor3_op_ifi(fld2,fld9,lx9,ly9,lz9,ixm92
     $      ,iytm92,iztm92,lx2,ly2,lz2,ifx,ify,ifz,tnsr_wk,lx1)
      elseif (idir.eq.-1) then
        call tensor3_op_ifi(fld9,fld2,lx2,ly2,lz2,ixm29
     $      ,iytm29,iztm29,lx9,ly9,lz9,ifx,ify,ifz,tnsr_wk,lx1)
      else
        if (nid.eq.0) write(6,*) 'Unknown idir in map92.'
        call exitt  
      endif  

      return
      end subroutine
!----------------------------------------------------------------------

      subroutine geodat91

C     Routine to generate elemental geometric matrices on Mesh 1
C     (Gauss-Legendre Lobatto mesh), based on geometry on Mesh 9

      implicit none  

      INCLUDE 'SIZE'
      INCLUDE 'GEOM'
      INCLUDE 'INPUT'
      INCLUDE 'MASS'
      INCLUDE 'TSTEP'
      INCLUDE 'WZ'

C     Note: Subroutines GLMAPM1, GEODAT1, AREA2, SETWGTR and AREA3 
C           share the same array structure in Scratch Common /SCRNS/.

!     These variables are calculated in glmapm9      
      REAL XRM1,YRM1,XSM1,YSM1,XTM1,YTM1,ZRM1,ZSM1,ZTM1,WJ
      COMMON /SCRNS/ XRM1(LX1,LY1,LZ1,LELT)
     $ ,             YRM1(LX1,LY1,LZ1,LELT)
     $ ,             XSM1(LX1,LY1,LZ1,LELT)
     $ ,             YSM1(LX1,LY1,LZ1,LELT)
     $ ,             XTM1(LX1,LY1,LZ1,LELT)
     $ ,             YTM1(LX1,LY1,LZ1,LELT)
     $ ,             ZRM1(LX1,LY1,LZ1,LELT)
      COMMON /CTMP1/ ZSM1(LX1,LY1,LZ1,LELT)
     $ ,             ZTM1(LX1,LY1,LZ1,LELT)
     $ ,             WJ   (LX1,LY1,LZ1,LELT)

      integer iel,nxyz1,ntot1
      integer i,j

      nxyz1 = lx1*ly1*lz1
      ntot1 = nxyz1*nelt

!      IF (IFGMSH3 .AND. ISTEP.EQ.0)
!     $   CALL XYZRST (XRM1,YRM1,ZRM1,XSM1,YSM1,ZSM1,XTM1,YTM1,ZTM1,
!     $                IFAXIS)

      if (.not.ifaxis) then
        call invers2 (wj,jacm1,ntot1)
      else
        do iel=1,nelt
          if (ifrzer(iel)) then
            do j=1,ly1
            do i=1,lx1
              if (j.gt.1) then
                 wj(i,j,1,iel) = ym1(i,j,1,iel)/
     $                          (jacm1(i,j,1,iel)*(1.+zam1(j)))
              else
                 wj(i,j,1,iel) = ysm1(i,j,1,iel)/jacm1(i,j,1,iel)
              endif
            enddo
            enddo
          else
            call invcol3 (wj(1,1,1,iel),ym1(1,1,1,iel),
     $                    jacm1(1,1,1,iel),nxyz1)
          endif
        enddo
      endif

!     Compute geometric factors for integrated del-squared operator.

      if (ldim.eq.2) then
        call vdot2 (g1m1,rxm1,rym1,rxm1,rym1,ntot1)
        call vdot2 (g2m1,sxm1,sym1,sxm1,sym1,ntot1)
        call vdot2 (g4m1,rxm1,rym1,sxm1,sym1,ntot1)
        call col2  (g1m1,wj,ntot1)
        call col2  (g2m1,wj,ntot1)
        call col2  (g4m1,wj,ntot1)
        call rzero (g3m1,ntot1)
        call rzero (g5m1,ntot1)
        call rzero (g6m1,ntot1)
      else
        call vdot3 (g1m1,rxm1,rym1,rzm1,rxm1,rym1,rzm1,ntot1)
        call vdot3 (g2m1,sxm1,sym1,szm1,sxm1,sym1,szm1,ntot1)
        call vdot3 (g3m1,txm1,tym1,tzm1,txm1,tym1,tzm1,ntot1)
        call vdot3 (g4m1,rxm1,rym1,rzm1,sxm1,sym1,szm1,ntot1)
        call vdot3 (g5m1,rxm1,rym1,rzm1,txm1,tym1,tzm1,ntot1)
        call vdot3 (g6m1,sxm1,sym1,szm1,txm1,tym1,tzm1,ntot1)
        call col2  (g1m1,wj,ntot1)
        call col2  (g2m1,wj,ntot1)
        call col2  (g3m1,wj,ntot1)
        call col2  (g4m1,wj,ntot1)
        call col2  (g5m1,wj,ntot1)
        call col2  (g6m1,wj,ntot1)
      endif

!     Multiply the geometric factors GiM1,i=1,5 with the
!     weights on mesh M1.

      do iel=1,nelt
        if (ifaxis) call setaxw1 ( ifrzer(iel) )
          call col2 (g1m1(1,1,1,iel),w3m1,nxyz1)
          call col2 (g2m1(1,1,1,iel),w3m1,nxyz1)
          call col2 (g4m1(1,1,1,iel),w3m1,nxyz1)
        if (ldim.eq.3) then
          call col2 (g3m1(1,1,1,iel),w3m1,nxyz1)
          call col2 (g5m1(1,1,1,iel),w3m1,nxyz1)
          call col2 (g6m1(1,1,1,iel),w3m1,nxyz1)
        endif
      enddo

!     Compute the mass matrix on mesh M1.

      do iel=1,nelt
        if (ifaxis) call setaxw1 ( ifrzer(iel) )
        call col3 (bm1  (1,1,1,iel),jacm1(1,1,1,iel),w3m1,nxyz1)
        if (ifaxis) then 
          call col3(baxm1(1,1,1,iel),jacm1(1,1,1,iel),w3m1,nxyz1)
          if (ifrzer(iel)) then
            do j=1,ly1
              if (j.gt.1) then
                 do i=1,lx1
                   bm1(i,j,1,iel) = bm1(i,j,1,iel)*ym1(i,j,1,iel)
     $                                            /(1.+zam1(j))
                   baxm1(i,j,1,iel)=baxm1(i,j,1,iel)/(1.+zam1(j))
                 enddo  
              else
                 do i=1,lx1
                   bm1(i,j,1,iel) = bm1(i,j,1,iel)*ysm1(i,j,1,iel)
                   baxm1(i,j,1,iel)=baxm1(i,j,1,iel)
                 enddo  
              endif
            enddo
          else
            call col2 (bm1(1,1,1,iel),ym1(1,1,1,iel),nxyz1)
          endif
        endif
      enddo

      if(ifaxis) then
        do iel=1,nelt
          if(ifrzer(iel)) then
            do j=1,ly1
            do i=1,lx1
              if(j.eq.1) then
                 yinvm1(i,j,1,iel)=1.0/ysm1(i,j,1,iel)
              else
                 yinvm1(i,j,1,iel)=1.0/ym1 (i,j,1,iel)
              endif
            enddo 
            enddo 
          else
            call invers2(yinvm1(1,1,1,iel),ym1(1,1,1,iel),nxyz1)
          endif
        enddo
      else
        call cfill(yinvm1,1.0,nxyz1*nelt)
      endif

!     Compute normals, tangents, and areas on elemental surfaces
!     Relies of dx/dr etc. on Mesh 1.
!     We have calculated those based on low-order Mesh 9
      call setarea

      return
      end
!---------------------------------------------------------------------- 

      subroutine updcoorM9

!     Subroutine to update geometry for moving boundary problems

      implicit none

      include 'SIZE'
      include 'TSTEP'
      include 'INPUT'

      integer nel

      ifield = 0
      nel    = nelfld(ifield)

!     Update collocation points coordinates

      call updxyzM9 (nel)

!     Shift lagged mesh velocity

      if (.not.ifrich) call lagmshvM9 (nel)

      return
      end
!-----------------------------------------------------------------------

      subroutine updxyzM9 (nel)

      implicit none

      include 'SIZE'
      include 'TSTEP'
      include 'MVGEOM'
      include 'GEOM'
      include 'INPUT'

      include 'SUBGEOM'

      integer nel

      real ux,uy,uz
      common /scrsf/ ux(lx9,ly9,lz9,lelt)
     $             , uy(lx9,ly9,lz9,lelt)
     $             , uz(lx9,ly9,lz9,lelt)
      real abm(3)

      integer i,ilag,ntot9

      ntot9 = lx9*ly9*lz9*nel

      do i=1,nbd
        abm(i) = dt*abmsh(i)
      enddo  

      if (istep.eq.0) then
        call copy (ux,wx9,ntot9)
        call copy (uy,wy9,ntot9)
        if (ldim.eq.3) call copy (uz,wz9,ntot9)
      else
        if (ifrich) then
          call cmult2(ux,wx,dt,ntot9)
          call cmult2(uy,wy,dt,ntot9)
          if (ldim.eq.3) call cmult2(uz,wz,dt,ntot9)
        else
          call cmult2 (ux,wx,abm(1),ntot9)
          call cmult2 (uy,wy,abm(1),ntot9)
          if (ldim.eq.3) call cmult2 (uz,wz,abm(1),ntot9)
          do ilag=2,nbd
            call add2s2 (ux,wxlag9(1,1,1,1,ilag-1),abm(ilag),ntot9)
            call add2s2 (uy,wylag9(1,1,1,1,ilag-1),abm(ilag),ntot9)
            if (ldim.eq.3)
     $      call add2s2 (uz,wzlag9(1,1,1,1,ilag-1),abm(ilag),ntot9)
          enddo   
        endif
      endif

      call add2 (xm9,ux,ntot9)
      call add2 (ym9,uy,ntot9)
      if (ldim.eq.3) call add2 (zm9,uz,ntot9)

      return
      end
!-----------------------------------------------------------------------
      subroutine lagmshvM9 (nel)

!     Keep old mesh velocity

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'MVGEOM'
      include 'TSTEP'
      include 'SUBGEOM'
      
      integer nel,il,ntot9

      ntot9 = lx9*ly9*lz9*nel

      do il=nbdinp-1,2,-1
         call copy (wxlag9(1,1,1,1,il),wxlag9(1,1,1,1,il-1),ntot9)
         call copy (wylag9(1,1,1,1,il),wylag9(1,1,1,1,il-1),ntot9)
         if (ldim.eq.3)
     $   call copy (wzlag9(1,1,1,1,il),wzlag9(1,1,1,1,il-1),ntot9)
      enddo

      call copy (wxlag9(1,1,1,1,1),wx9,ntot9)
      call copy (wylag9(1,1,1,1,1),wy9,ntot9)
      if (ldim.eq.3) call copy (wzlag9(1,1,1,1,1),wz9,ntot9)

      return
      end
c-----------------------------------------------------------------------






