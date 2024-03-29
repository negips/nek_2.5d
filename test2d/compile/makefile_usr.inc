
SRC_INC     := ./src_inc
SRC_USR     := ./src_usr
SRC_FRAME   := ./frame
SRC_ARN     := ./src_arn

SRC_F3D     := ./src_f3d
SRC_FS      := ./src_fs
SRC_CHP     := ./checkpoint

$(OBJDIR)/frame_usr.o           :$(SRC_FRAME)/frame_usr.f;              $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/frame.o               :$(SRC_FRAME)/frame.f;                  $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/mntrlog_block.o       :$(SRC_FRAME)/mntrlog_block.f;          $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/mntrlog.o             :$(SRC_FRAME)/mntrlog.f;                $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/mntrtmr_block.o       :$(SRC_FRAME)/mntrtmr_block.f;          $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/mntrtmr.o             :$(SRC_FRAME)/mntrtmr.f;                $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/rprm_block.o          :$(SRC_FRAME)/rprm_block.f;             $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/rprm.o                :$(SRC_FRAME)/rprm.f;                   $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/io_tools.o            :$(SRC_FRAME)/io_tools.f;               $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/io_tools_block.o      :$(SRC_FRAME)/io_tools_block.f;         $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/math_tools.o          :$(SRC_FRAME)/math_tools.f;             $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/solve_f3d.o           :$(SRC_F3D)/solve_f3d.f ;               $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/solvep_f3d.o          :$(SRC_F3D)/solvep_f3d.f;               $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/weak_lapl_cyl.o       :$(SRC_F3D)/weak_lapl_cyl.f;            $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/solver_f3d.o          :$(SRC_F3D)/solver_f3d.f;               $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/frame_f3d.o           :$(SRC_F3D)/frame_f3d.f;                $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/tstepper.o            :$(SRC_ARN)/tstepper.f;                 $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/conht_tools.o         :$(SRC_ARN)/conht_tools.f;              $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/arn_arp.o             :$(SRC_ARN)/arn_arp.f;                  $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/arn_arp_io.o          :$(SRC_ARN)/arn_arp_io.f;               $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/chkpoint.o            :$(SRC_CHP)/chkpoint.f;                 $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/chkpt_mstp.o          :$(SRC_CHP)/chkpt_mstp.f;               $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/arn_cmath.o           :$(SRC_ARN)/arn_cmath.f;                $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@

$(OBJDIR)/fs_ale.o              :$(SRC_FS)/fs_ale.f;                    $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/fs_sort.o             :$(SRC_FS)/fs_sort.f;                   $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/fs_sort_mp.o          :$(SRC_FS)/fs_sort_mp.f;                $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/fs_global_surf.o      :$(SRC_FS)/fs_global_surf.f;            $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/fs_global_gmres.o     :$(SRC_FS)/fs_global_gmres.f;           $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/frame_fs.o            :$(SRC_FS)/frame_fs.f;                  $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/localfilter.o         :$(SRC_USR)/localfilter.f;              $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@

$(OBJDIR)/generic_subs.o        :$(SRC_USR)/generic_subs.f;             $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/kopriva.o             :$(SRC_USR)/kopriva.f;                  $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/subparametric_geom.o  :$(SRC_USR)/subparametric_geom.f;       $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@


