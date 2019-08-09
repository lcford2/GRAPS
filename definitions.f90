
!   Declaration of prototype variables

!   Declaration of System Blocks and Demand Descriptors

 MODULE  Definitions

	TYPE User

		Character*40 Name
		Integer   ID,nchild,nparent,user_type,nres_level,nlags
		Integer, allocatable ::  parent_id(:),child_id(:)
		Integer, allocatable ::  parent_type(:),child_type(:)
		Integer, allocatable :: ffraction(:)

		double precision minimum_release, maximum_release
		double precision tariff,penalty, failure_prob,con_res_vol,penalty_compen
		double precision, allocatable :: demand_fract(:), res_compensation(:),restr_fract(:)

		double precision max_discharge,installed_capacity,generator_efficiency,storage_energy_coeff(2)
		double precision unit_rate_energy
                double precision, allocatable :: tail_elevation(:)



   END TYPE User

	TYPE Reservoir

		Character*40 Name
		Integer   ID
		Integer	 nspillway
		Integer  number_outlets
		Integer  nres_level
		Integer  nchild, nparent

		Integer, allocatable ::  parent_id(:),child_id(:),spill_type(:)
		Integer, allocatable ::  parent_type(:),child_type(:)


		double precision latitude, longitude
		double precision elevation_max, elevation_min, storage_max, storage_min
		double precision elevation_storage_coeff(3), storage_area_coeff(2)

		double precision, allocatable :: elevation_rvalve(:), area_rvalve(:),loss_coeff_max(:), loss_coeff_min(:) 
		double precision current_storage, target_storage, storage_prob
		double precision, allocatable ::  rule_curve(:),evaporation_rate(:)
		double precision, allocatable :: tar_restr_prob(:),restr_fraction(:),crest_level(:),discharge_max(:)
		

	END TYPE Reservoir

   TYPE Flow_join_node

		Character*40 Name
		Integer  nuser
		Integer   ID,nchild,nparent
		Integer, allocatable ::  parent_id(:),child_id(:)
		Integer, allocatable ::  parent_type(:),child_type(:)

  END TYPE Flow_join_node


TYPE Sink

		Character*40 Name
		Integer   ID,nparent

		Integer, allocatable ::  parent_id(:)
		Integer, allocatable ::  parent_type(:)

		double precision max_storage


END TYPE Sink


TYPE Watershed

	Character*40 Name
	Integer   ID,nchild

	Integer, allocatable ::  child_id(:)
	Integer, allocatable ::  child_type(:)

	double precision drainage_area
	double precision, allocatable :: natural_inflows(:,:)

END TYPE Watershed

TYPE Interbasin

	Character*40 Name
	Integer   ID,nchild

	Integer, allocatable ::  child_id(:)
	Integer, allocatable ::  child_type(:)

	double precision drainage_area
	double precision, allocatable :: average_flow(:)


END TYPE Interbasin


!   Declarion of Flow Connectors

TYPE direct_inflows

	Character*40 Name
	Integer start_id, end_id
	Integer start_type, end_type

	double precision minimum_discharge, maximum_discharge,loss_factor

END TYPE direct_inflows

TYPE Diversion

	Character*40 Name
	Integer start_id, end_id
	Integer start_type, end_type

	double precision minimum_discharge, maximum_discharge,loss_factor

END TYPE Diversion

TYPE passage_flow

	Character*40 Name
	Integer start_id, end_id
	Integer start_type, end_type

	double precision minimum_discharge, maximum_discharge,loss_factor

END TYPE passage_flow


TYPE return_inflows

	Character*40 Name
	Integer start_id, end_id
	Integer start_type, end_type

	double precision minimum_discharge, maximum_discharge,loss_factor

END TYPE return_inflows


TYPE demand_release

	Character*40 Name
	Integer start_id, end_id
	Integer start_type, end_type

	double precision minimum_discharge, maximum_discharge,loss_factor

END TYPE demand_release



TYPE Interbasin_flow

	Character*40 Name
	Integer start_id, end_id
	Integer start_type, end_type

	double precision minimum_discharge, maximum_discharge,loss_factor

END TYPE Interbasin_flow

TYPE Natural_flow

	Character*40 Name
	Integer start_id, end_id
	Integer start_type, end_type

	double precision minimum_discharge, maximum_discharge,loss_factor

END TYPE Natural_flow

TYPE Spill_flow

	Character*40 Name
	Integer start_id, end_id
	Integer start_type, end_type

	double precision minimum_discharge, maximum_discharge,loss_factor

END TYPE Spill_flow

TYPE ordered_network

	Integer order_id, order_type

END TYPE ordered_network

TYPE flow_definitions

	doubleprecision, allocatable :: controlled_flows(:)
	doubleprecision, allocatable :: uncontrolled_flows(:,:)

END TYPE flow_definitions

END MODULE Definitions
