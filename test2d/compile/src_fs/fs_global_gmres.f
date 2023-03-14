!---------------------------------------------------------------------- 
      subroutine fs_gmres(r1,r2,maxiter)

!     using right-preconditioned GMRES iteration.
!     to solve \rho*B*u = f

!     We solve:   \rho*B*(M^-1)(Mu) = f,
!              => \rho*B*(M^-1)v    = f,
!              => u = (M^-1)v

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'TSTEP'
      
      include 'FS_ALE'
      include 'FS_GMRES'


      integer ltfs
      parameter (ltfs = lxfs*lyfs)

      logical          ifprint
      common  /cprint/ ifprint

      real             r1    (ltfs)
      real             r2    (ltfs)
      real             wk    (ltfs)

      real theta        ! Approximated eigenvalue

      real alpha, l, temp
      integer j,m

      logical iflag
      save    iflag
      data    iflag /.false./
      real    norm_fac
      save    norm_fac

      real*8 etime1,dnekclock

      integer nt1,nt2,nt3
      real glsc2,glsc3,vlsc2,vlsc3,op_glsc2_wt
      integer iconv,intype
      real tolpss,div0
      integer i,k,iter
      real etime2,etime_p,ratio,rnorm

      integer maxiter

      logical ifwgt           ! If weighted Gmres       
      logical ifprec

      integer ngs             ! No of Gram-Schmidt
      character*132 str

      real D1(ltfs),D2(ltfs),D3(ltfs)

      integer i1,i2,i3


      nt1      = lxfs*lyfs
      nt2      = nt1*2

      i1       = 1
      i2       = 1 + nt1
!      i3       = 1 + 2*nt1

      ifprint  = .true.

      ifwgt    = .false.           ! Weighted Orthogonalization
      ifprec   = .false.           ! Use preconditioner
      ngs      = 1

!     I use this for weights
      if (ifwgt) then
!        call opcopy(wghts(i1),wghts(i2),wghts(i3),bm1,bm1,bm1)
        call copy(wghts(i1),bm_fs,nt1)
        call copy(wghts(i2),bm_fs,nt1)
      else
        call rone(wghts,nt2)
      endif  
      alpha = sqrt(glsc2(wghts,wghts,nt2))
      norm_fac = 1.0/alpha

C     Set up diag preconditioner.
      if (ifprec) then
        call rone(D1,nt1)  
        call rone(D2,nt1)  
        call rone(D3,nt1)  
      else
        call rone(D1,nt1)  
        call rone(D2,nt1)  
        call rone(D3,nt1)  
      endif  

      etime1 = dnekclock()
      etime_p = 0.
      iter  = 0
      m = min(maxiter,lgmit)

      tolpss = 1.0e-12

      iconv = 0
      call rzero   (x_gmres,ltgm)

      do while(iconv.eq.0.and.iter.lt.maxiter)

         if (iter.eq.0) then
!            call opcopy(r_gmres(i1),r_gmres(i2),r_gmres(i3),
!     $                  r1,r2,r3)
           call copy(r_gmres(i1),r1,nt1)
           call copy(r_gmres(i2),r2,nt1)
         else
!           update residual
!            call opcopy(r_gmres(i1),r_gmres(i2),r_gmres(i3),
!     $                  r1,r2,r3)
           call copy(r_gmres(i1),r1,nt1)
           call copy(r_gmres(i2),r2,nt1)

!          w = A*x
           call fs_lagrangian_deriv(w_gmres(i1),w_gmres(i2),
     $                              x_gmres(i1),x_gmres(i2),wk)

!          r = r - w
           call add2s2(r_gmres,w_gmres,-1.,nt2)

         endif

         if (ifwgt) then
!          Weighted inner product                                
!                     ______
!          gamma  = \/(Br,r) 
           gamma_gmres(1) = sqrt(glsc3(r_gmres,r_gmres,wghts,nt2))
         else    
!          Un-weighted inner product
!                     ______
!          gamma  = \/(r,r)                          
           gamma_gmres(1) = sqrt(glsc2(r_gmres,r_gmres,nt2))
         endif   
                                                           
         if(iter.eq.0) then
           div0 = gamma_gmres(1)*norm_fac
           if (param(21).lt.0) tolpss=abs(param(21))*div0
         endif

!        check for lucky convergence
         rnorm = 0.
         if (gamma_gmres(1) .eq. 0.) then
           iconv = 1 
           exit   ! exit outer loop
         endif  
!goto 9000
         temp = 1./gamma_gmres(1)
         call cmult2(v_gmres(1,1),r_gmres,temp,nt2)! v  = r / gamma

         do j=1,m
            iter = iter+1

            call copy(w_gmres,v_gmres(1,j),nt2) ! w  = v_j

            etime2 = dnekclock()
!           z = (M^-1)w      
            if (ifprec) then
              call col3(z_gmres(i1,j),w_gmres(i1),D1,nt1)
              call col3(z_gmres(i2,j),w_gmres(i2),D2,nt1)
            else
              call copy(z_gmres(1,j),w_gmres,nt2)
            endif

            etime_p = etime_p + dnekclock()-etime2

!           r = (M^-1)w 
            call copy(r_gmres,z_gmres(1,j),nt2)

!           w = A*(M^-1)w
            call fs_lagrangian_deriv(w_gmres(i1),w_gmres(i2),
     $            r_gmres(i1),r_gmres(i2),wk)

!           Gram-Schmidt:
            call ortho_subspace(w_gmres,nt2,h_gmres(1,j),v_gmres,
     $            ltgm,j,wghts,ifwgt,ngs,wk1,wk2,.true.)

!           Apply Givens rotations to new column
            do i=1,j-1
              temp = h_gmres(i,j)                   
              h_gmres(i  ,j)=  c_gmres(i)*temp 
     $                       + s_gmres(i)*h_gmres(i+1,j)  
              h_gmres(i+1,j)= -s_gmres(i)*temp 
     $                       + c_gmres(i)*h_gmres(i+1,j)
            enddo

            if (ifwgt) then
!                        ______
!             alpha =  \/(Bw,w) 
              alpha = sqrt(glsc3(w_gmres,w_gmres,wghts,nt2))
            else
!                        ______
!             alpha =  \/(w,w) 
              alpha = sqrt(glsc2(w_gmres,w_gmres,nt2))        
            endif
            rnorm = 0.

            if (alpha.eq.0.) then
              iconv = 1
              exit      ! exit inner loop
            endif
!            write(6,*) 'h_gmres(:,j)', (h_gmres(k,j),k=1,j)
            l = sqrt(h_gmres(j,j)*h_gmres(j,j)+alpha*alpha)
            temp = 1./l
            c_gmres(j) = h_gmres(j,j) * temp
            s_gmres(j) = alpha  * temp
            h_gmres(j,j) = l
            gamma_gmres(j+1) = -s_gmres(j) * gamma_gmres(j)
            gamma_gmres(j)   =  c_gmres(j) * gamma_gmres(j)
           
            rnorm = abs(gamma_gmres(j+1))*norm_fac
            ratio = rnorm/div0
            if (ifprint.and.nio.eq.0) 
     $         write (6,66) iter,tolpss,rnorm,div0,ratio,istep
   66       format(i5,1p4e12.5,i8,' FS Saddle ')

            if (rnorm .lt. tolpss) then
              iconv = 1
              exit ! exit inner loop
            endif   

            if (j.eq.m) goto 1000 !not converged, restart

            temp = 1./alpha
            call cmult2(v_gmres(1,j+1),w_gmres,temp,nt2) ! v = w / alpha
         enddo

 1000    continue
         !back substitution
         !     -1
         !c = H   gamma
         do k=j,1,-1
           temp = gamma_gmres(k)
           do i=j,k+1,-1
             temp = temp - h_gmres(k,i)*c_gmres(i)
           enddo
           c_gmres(k) = temp/h_gmres(k,k)
         enddo
!        sum up Arnoldi vectors
!        x_gmres = (M^-1)*(V*c)
!     => x_gmres = (M^-1 * V)*c = Z*c
         do i=1,j
!          x = x + Z*c
           call add2s2(x_gmres,z_gmres(1,i),c_gmres(i),nt2) 
         enddo
      enddo       ! outer loop
! 9000 continue

      call copy(r1,x_gmres(i1),nt1)
      call copy(r2,x_gmres(i2),nt1)

      etime1 = dnekclock()-etime1
      if (nio.eq.0) write(6,9999) istep,' FS Saddle gmres  ', 
     &                            iter,rnorm,div0,tolpss,etime_p,etime1

 9999 format(i11,a,I6,1p5e13.4)

      return
      end

!-----------------------------------------------------------------------

      subroutine fs_lagrangian_deriv(DLr1,DLr2,r1,r2,wk)

      implicit none

      include 'SIZE'
      include 'FS_ALE'

      real r1(lxfs,lyfs),r2(lxfs,lyfs)
      real DLr1(lxfs,lyfs),DLr2(lxfs,lyfs)
      real wk(lxfs*lyfs)

      integer i,j,k,n
      
      n = lxfs*lyfs

!     Interior lagrange multipliers are zero   
      do j=2,lyfs
      do i=1,lxfs
        r2(i,j) = 0.0
      enddo
      enddo  
     
!     Transposed operator here
      call copy(DLr1,r2,n)
      call col2(DLr1,w2_fs,n)
      call tensory_op(wk,DLr1,lxfs,lyfs,1,dy_fs,lyfs)

      call copy(DLr1,r1,n)
      call col2(DLr1,bm_fs,n)
      call add2(DLr1,wk,n)

!     Regular gradient operator here      
      call tensory_op(DLr2,r1,lxfs,lyfs,1,dyt_fs,lyfs)
      call col2(DLr2,w2_fs,n)
!     Interior lagrange multipliers are zero   
      do j=2,lyfs
      do i=1,lxfs
        DLr2(i,j) = 0.0
      enddo
      enddo  

      return
      end subroutine
!---------------------------------------------------------------------- 

      subroutine fs_lagrangian_deriv1D(DLr1,DLr2,r1,r2,wk)

      implicit none

      include 'SIZE'
      include 'FS_ALE'

      real r1(lxfs,1),r2(lxfs,1)
      real DLr1(lxfs,lyfs),DLr2(lxfs,lyfs)
      real wk(lxfs*lyfs)

      integer i,j,k,n
      
      n = lxfs*1

!     Interior lagrange multipliers are zero   
      do j=1,1
      do i=2,lxfs
        r2(i,j) = 0.0
      enddo
      enddo  

     
!     Transposed operator here
      call copy(DLr1,r2,n)
      call col2(DLr1,wx_fs,n)
      call tensorx_op(wk,DLr1,lxfs,1,1,dxt_fs,lxfs)

      call copy(DLr1,r1,n)
      call col2(DLr1,bm_fs,n)
!      call col2(DLr1,jac_fs,n)
      call add2(DLr1,wk,n)

!     Regular gradient operator here      
      call tensorx_op(DLr2,r1,lxfs,1,1,dx_fs,lxfs)
      call col2(DLr2,bm_fs,n)
!     Interior lagrange multipliers are zero   
      do j=1,1 !lyfs
      do i=2,lxfs
        DLr2(i,j) = 0.0
      enddo
      enddo  

!!     Dirichlet      
!      do j=1,1
!      do i=1,lxfs
!        DLr1(i,j) = 0.0
!      enddo
!      enddo  


      return
      end subroutine
!---------------------------------------------------------------------- 






