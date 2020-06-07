MODULE my_variables

USE DEFINITIONS

TYPE (User), allocatable			 ::	my_user(:)
TYPE (Reservoir), allocatable		 ::	my_reservoir(:)
TYPE (Flow_join_node), allocatable	 ::	my_node(:)
TYPE (Sink), allocatable			 ::	my_sink(:)
TYPE (Watershed), allocatable		 ::	my_watershed(:)
TYPE (direct_inflows), allocatable	 ::	my_dir_inflows(:)
TYPE (return_inflows), allocatable	 ::	my_ret_inflows(:)
TYPE (diversion), allocatable		 ::	my_diversions(:)
TYPE (Spill_flow), allocatable		 ::	my_spill_flow(:)
TYPE (Natural_flow), allocatable	 ::	my_natural_flow(:)
TYPE (Interbasin_flow), allocatable  ::	my_interbasin_flow(:)
TYPE (demand_release), allocatable	 ::	my_demand_release(:)
TYPE (Interbasin), allocatable		 ::	my_interbasin(:)
TYPE (ordered_network), allocatable  :: my_network_order(:),searched_vertices(:),parallel_track(:)
TYPE (flow_definitions), allocatable :: my_flow_set(:)
character*40 type_details
character*25, objective_function

Integer itype,nres,nuser,nfnode,nsink,nwatershed,ndir_inflows,nret_inflows,ndiversion,nspill_flow,nnatural_flow,ninterbasin_flow,ndemand_release
Integer ntime,nensem,ncons,iblock,ntotal_vertices,isimul_block,nsimul_block,iflow_set,nres_level,nperiods

Integer, allocatable :: isimul_status(:)

double precision, allocatable ::	temp_decision_var(:),cons_global(:),value_net(:)
double precision, allocatable ::    spill_values(:)


END MODULE my_variables