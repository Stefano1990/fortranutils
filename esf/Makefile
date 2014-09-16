# LESter

PROGRAM = LESter_master_output.prog 

SOURCE = main.F \
	setupcomm.F set_offset.F read_bl_size.F read_bl_info.F \
	extr_string.F set_pointer.F read_xyz.F mg_xyz.F xyz_extrap0.F \
	set_const.F read_3.F file_3.F \
	comp_metric.F metric_block2.F init_array.F swap_xyz.F \
	set_loop_bound.F metric_swap1.F metric_swap2.F file_2.F subgrid.F \
	subgrid_nosgs.F subgrid_wale.F subgrid_smagorinsky.F \
        subgrid_sigma.F \
        set_virt_bc.F set_virt_bc_block.F swap_sc0.F \
	set_uv_int.F set_u_int_block.F set_w_int.F \
	set_w_int_block.F coeff_p.F coeff_p_block.F sour_p.F sour_p_block.F \
	solvers.F lmg_solve.F norms.F norms_block.F solversc.F smooth.F \
	restriction.F prolongation.F sporlx.F rest_block.F prol_block.F \
	uvw_adv.F uvw_adv_block.F cflest.F cflest_block.F extrap_p.F \
	extrap_p_block.F p_ref_subtract.F global_sum.F file_2b.F \
	read_2b.F read_twall.F file_twall.F errhan_dummy.F \
	set_sc1var.F set_sc1var_block.F \
	conv_exitbc1_block.F conv_exitbc2_block.F statistics.F \
	conv_exitbc_scale.F mass_flux.F \
	step.F conv_exitbc.F fixed_inlet.F initia.F \
	set_gamsc1_bc.F set_sc1_bc.F set_sc1.F set_sc1_block.F \
	set_interfaces.F sporlx2.F usr_sum.F \
	swap_rhos.F swap_rhos_block.F search_lib.F search_lib_block.F \
	mass_change.F mass_change_block.F defaults.F \
        dataread.F field_out.F assignfile.F slipbc.F

OBJEXT  = .o
OBJ = $(SOURCE:.F=.o)

#ARCHF = -DSGI -DMPI1 -O3 -xSSE2 -r8  
# Stephan: replaced -xSSE2 with -xSSE3
ARCHF = -DSGI -DMPI1 -O3 -xSSE3 -r8
FC=mpiifort


# Other flags for the preprocessor:
# -DMPI1   : Use only MPI-1.X standard routines, not MPI-2
# -DSHMEM  : Use SHMEM in swap_sc (add -lsma to LDFLAGS2)
# -DDEBUG  : Compile a version suitable for debugging (i.e.
#            one using ordered global sums).
#
# Note - you can have DEBUG and one of the other two, but
# not SHMEM and MPI1 at the same time. Also note that you
# don't need to turn on MPI1 on platforms where it is the
# default (see relevant ARCHF above). Of course, you can
# only turn on SHMEM if it's available on your platform
# (Cray T3E and SGI Origin should have this).

#OTHERF = -DMPI1
#OTHERF = -DSHMEM
#OTHERF = -DDEBUG

# You should not need to edit anything below this line
#------------------------------------------------------------


FFLAGS = $(ARCHF) $(OTHERF)


.SUFFIXES:
.SUFFIXES: $(OBJEXT) .F


$(PROGRAM) : $(OBJ)
	$(FC) $(LDFLAGS) $(OBJ) -o $@ $(LDFLAGS2)

.F$(OBJEXT):
	$(FC) -c $(FFLAGS) $< -o $@

clean:
	rm -f core \
	$(PROGRAM) \
	$(OBJ)

depend:
	./sfmakedepend -f Makedepend $(SOURCE)

include Makedepend
