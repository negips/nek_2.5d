!----------------------------------------------------------------------
!     Author: Prabal Negi
!     Description: Routines for Constrained motion of the interface.
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      subroutine fs_constr_meshmv(wx,wy,wz)

      implicit none

      include 'SIZE'
      include 'FS_ALE'
!      include 'SOLN'

      real wx(lx1,ly1,lz1,lelv),wy(lx1,ly1,lz1,lelv)
      real wz(lx1,ly1,lz1,lelv)

      integer icalld
      save icalld
      data icalld /0/

      integer i

!     Generate the global basis anyway      
      if (icalld.eq.0) then
        call fs_global_basis
        icalld = icalld+1
      endif

!     Zero out tangential component of mesh velocity
!      call fs_mvmeshn2(wx,wy,wz)


      if (fs_ifgsm) then

        call fs_gllo_xyz

        call fs_gllo_flds(wx,wy,wz)
        call fs_intp_setup
        call fs_get_localpts          ! SEM -> Global

!       Filtering here requires BOYD transformation to be active.
!       Otherwise the boundary points move.
!       Alternately, one could correct boundary points.      
!!       Filter normal velocities
!        if (fs_iffil) then
!          do i=1,ndim
!            call fs_glfilter(fld_fs(1,1,i))
!          enddo  
!        endif
        call fs_get_globalpts         ! Global -> SEM
        call fs_restore_int(wx,wy,wz,'Norm')
      endif       ! fs_ifgsm 

!     Correction for tangential movement to minimize
!     mesh deformation
      if (fs_iftc) then
        call fs_tang_corr(wx,wy,wz)
      endif

!     make sure fluid does not come off the wall      
      call fs_fixcorners(wx,wy,wz) 

!     Free the handles
      if (fs_ifgsm) then 
        call fgslib_findpts_free(intgh_fs)
        call fgslib_findpts_free(intlh_fs)

        call mntr_log(fs_id,fs_log,'Interface Smoothening Done')
      endif  

      return
      end subroutine fs_constr_meshmv
!----------------------------------------------------------------------

      subroutine fs_global_deriv()

      implicit none
        
      include 'SIZE'
      include 'FS_ALE'

     














      return
      end subroutine fs_global_deriv        
!---------------------------------------------------------------------- 
