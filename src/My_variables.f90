MODULE my_variables

USE DEFINITIONS

TYPE (User), allocatable			::	my_user(:)
TYPE (Reservoir), allocatable		::	my_reservoir(:)
TYPE (Flow_join_node), allocatable	::	my_node(:)
TYPE (Sink), allocatable			::	my_sink(:)
TYPE (Watershed), allocatable		::	my_watershed(:)
TYPE (direct_inflows), allocatable	::	my_dir_inflows(:)
TYPE (return_inflows), allocatable	::	my_ret_inflows(:)
TYPE (diversion), allocatable		::	my_diversions(:)
TYPE (Spill_flow), allocatable		::	my_spill_flow(:)
TYPE (Natural_flow), allocatable	::	my_natural_flow(:)
TYPE (Interbasin_flow), allocatable ::	my_interbasin_flow(:)
TYPE (demand_release), allocatable	::	my_demand_release(:)
TYPE (Interbasin), allocatable		::	my_interbasin(:)
TYPE (ordered_network), allocatable ::  my_network_order(:),searched_vertices(:),inproc_vertices(:),parallel_track(:)
TYPE (flow_definitions), allocatable ::  my_flow_set(:)
character*40 type_details
! character*500 input_path, output_path

integer itype,nres,nuser,nfnode,nsink,nwatershed,ndir_inflows,nret_inflows,ndiversion
integer nspill_flow,nnatural_flow,ninterbasin_flow,ndemand_release
integer ntime,nensem, ncons,iblock,ntotal_vertices,isimul_block,nsimul_block,iflow_set,nres_level,lukes_ncons,opt_count
double precision ben_net, constraint_tolerance
integer, allocatable :: isimul_status(:)

double precision, allocatable :: temp_decision_var(:),cons_global(:),value_net(:)

! variables for export to python
integer, allocatable :: id_output(:), cons_id(:), res_ids_for_spdef(:), user_id(:), cons_ignore(:)
double precision, allocatable :: value_output(:), spill_values(:), deficit_values(:)
double precision, allocatable :: cons_mag(:), min_rel(:), max_rel(:), constraints(:)
double precision, allocatable :: lukes_cons(:), hydro_benefit(:), all_release(:)

! inproc_vertices length, to treat it similar to a stack
integer :: inproc_top

END MODULE my_variables