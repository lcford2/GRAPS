!Program for multi-reservoir allocation
Program multireservoir
USE My_variables
implicit doubleprecision(a-h,o-z)
common /final_print/ifinal

external constr_res,expected_benefits,grobfd,grcnfd,ffsqp,simple_benefits,spill_objective,max_release

double precision, allocatable	::	decision_var_lb(:),decision_var_ub(:), decision_var(:),w(:),g(:),f(:)
Integer,allocatable :: iw(:)
double precision expect_bens, initial_bens


! Opening input files for reading
open(unit=10, file = 'input_data_files/input.dat',					ACTION = 'READ', STATUS = 'OLD', ERR = 100)
open(unit=11, file = 'input_data_files/watershed_details.dat',		ACTION = 'READ', STATUS = 'OLD')
open(unit=12, file = 'input_data_files/nflow_details.dat',			ACTION = 'READ', STATUS = 'OLD')
open(unit=13, file = 'input_data_files/reservoir_details.dat',		ACTION = 'READ', STATUS = 'OLD')
open(unit=14, file = 'input_data_files/user_details.dat',			ACTION = 'READ', STATUS = 'OLD')
open(unit=15, file = 'input_data_files/node_details.dat',			ACTION = 'READ', STATUS = 'OLD')
open(unit=16, file = 'input_data_files/dir_flow_details.dat',		ACTION = 'READ', STATUS = 'OLD')
open(unit=17, file = 'input_data_files/ret_flow_details.dat',		ACTION = 'READ', STATUS = 'OLD')
open(unit=18, file = 'input_data_files/diversions_details.dat',		ACTION = 'READ', STATUS = 'OLD')
open(unit=19, file = 'input_data_files/spill_flow_details.dat',		ACTION = 'READ', STATUS = 'OLD')
open(unit=20, file = 'input_data_files/ibasin_flow_details.dat',	ACTION = 'READ', STATUS = 'OLD')
open(unit=21, file = 'input_data_files/demand_flow_details.dat',	ACTION = 'READ', STATUS = 'OLD')
open(unit=22, file = 'input_data_files/sink_details.dat',			ACTION = 'READ', STATUS = 'OLD')
open(unit=23, file = 'input_data_files/interbasin_details.dat',		ACTION = 'READ', STATUS = 'OLD')
open(unit=24, file = 'input_data_files/runflag.dat',                ACTION = 'READ', STATUS = 'OLD')
open(unit=25, file=  'input_data_files/model_para.dat',				ACTION = 'READ', STATUS = 'OLD')
! Opening output files for writing
open(unit=31, file = 'output_files/storage.out')
open(unit=32, file = 'output_files/hydro.out')
open(unit=33, file = 'output_files/release.out')
open(unit=34, file = 'output_files/flow_sets.out')
open(unit=35, file = 'output_files/node_flow.out')
open(unit=36, file = 'output_files/spill.out')
open(unit=37, file = 'output_files/deficit.out')
open(unit=38, file = 'output_files/benefits.out')

!  reading input.dat 
read(10,*)ntime,nperiods,nensem
read(10,*)nwatershed,nnatural_flow,nres,nuser,nfnode,ndir_inflows,nret_inflows,ndiversion,nspill_flow,ninterbasin_flow,ndemand_release,nsink,ninterbasin

! checksum should equal the total number of system blocks and links. It will be the sum of the second line in `input.dat`
checksum = nwatershed+nnatural_flow+nres+nuser+nfnode+ndir_inflows+nret_inflows+ndiversion+nspill_flow+ninterbasin_flow+ndemand_release+nsink+ninterbasin

! ntime = Number of time steps
! nperiods = Number of times ntime will be repeated. Useful for multiyear runs at the monthly or season a scale.
!* 	nperiods is only used when reading in data, it is provided for user convenience to prevent extremely
!* 	long lines of data input (e.g. can stack inflow data in the file rather than providing it all on one line)
!* 	Example: 
!** 	run model for 4 years at monthly time steps. ntime = 12, nperiods = 4. 
!**		everywhere a monthly value is provided (et depth, inflow, rule curves, etc)
!** 	instead of providing 48 values on one line, provide 4 lines in succession
!**  	of 12 values per line. 
!* 	After input data is read, ntime will become ntime*nperiods and all output will be formatted as such
! nres = Number of reservoirs
! nuser = Number of water users in the entire system
! nfnode = Number of flow diversion/connection node
! nsink = Number of sink points or netowrks ends
! nwatershed - Number of watershed originating natural flows
! ndir_inflows - Number of Direct Inflows released from the reservoir
! nret_flows = Number of return flows (should not exceed nuser)
! ndiversion = Number of diversions towards environmental protection
! nspill_flow = Number of Spillways (should not exceed nres)
! nnatural_flow = Number of Natural flows from each watershed (should be equal to nwatershed)
! ninterbasin_flow = Number of Inter basin transfer points
! ndemand_release = Number of demand releases

! Allocate space for relevant data structures.

Allocate(my_user(nuser),my_reservoir(nres),my_node(nfnode),my_sink(nsink),my_watershed(nwatershed),my_interbasin(ninterbasin))

Allocate(my_dir_inflows(ndir_inflows),my_ret_inflows(nret_inflows),my_diversions(ndiversion))

Allocate(my_spill_flow(nspill_flow),my_natural_flow(nnatural_flow),my_interbasin_flow(ninterbasin_flow),&
    	 my_demand_release(ndemand_release))

! Reads the system details from individual files until the connectivity given in input.dat comes to an end.
icount = 0
icount_max = checksum

DO WHILE (icount<icount_max)
! call each subroutines and read input files 
  read(10,21)itype,type_details
	if(itype.eq.1)	call read_watershed_details(my_watershed,nwatershed,ntime,nensem,nperiods)
	if(itype.eq.2)	call read_nflow_details(my_natural_flow,nnatural_flow)
	if(itype.eq.3)	call read_reservoir_details(my_reservoir,nres,ntime,nensem,nperiods)
	if(itype.eq.4)	call read_user_details(my_user,nuser,ntime,nensem,nperiods)
  	if(itype.eq.5)	call read_node_details(my_node,nfnode)
	if(itype.eq.6)	call read_dir_inflows_details(my_dir_inflows,ndir_inflows)
	if(itype.eq.7)	call read_ret_inflows_details(my_ret_inflows,nret_inflows)
	if(itype.eq.8)	call read_diversions_details(my_diversions,ndiversion)
	if(itype.eq.9)	call read_spillflow_details(my_spill_flow,nspill_flow)
	if(itype.eq.10)	call read_interbasin_flow_details(my_interbasin_flow,ninterbasin_flow)
	if(itype.eq.11)	call read_demandrelease_details(my_demand_release,ndemand_release)
	if(itype.eq.12)	call read_sink_details(my_sink,nsink)
	if(itype.eq.13)	call read_interbasin_details(my_interbasin,ninterbasin,ntime,nperiods)

icount = icount +1 
END DO

! Close the input files
close(10)
close(11)
close(12)
close(13)
close(14)
close(15)
close(16)
close(17)
close(18)
close(19)
close(20)
close(21)
close(22)
close(23)

!   Convert the data structure into input parameters for the optimization routine.
ntime = ntime*nperiods
nparam = (nuser*ntime) ! each reservoir should be given a minimum of one user for downstream release and it should be defined in the user_details.dat
nsimul_block = nres + nfnode
ntotal_vertices = nwatershed + nres + nuser + nfnode + ninterbasin_flow + nsink 

Allocate(decision_var_lb(nparam), decision_var_ub(nparam), decision_var(nparam))

Allocate(spill_values(nparam))

Allocate(temp_decision_var(nparam),isimul_status(nsimul_block),my_network_order(nsimul_block))

Allocate(searched_vertices(ntotal_vertices))

! Initialize temporary decision variable
call array_ini(nparam,temp_decision_var,-10.0d0)

! Initialize spill values for optimization
call array_ini(nparam, spill_values, 0.0)

call solution_path()

! print *, (k1, my_network_order(k1)%order_type, my_network_order(k1)%order_id, k1=1, nsimul_block)

! Extract lower and upper bounds of decision variables.

call constr_decision_extract(decision_var_lb,decision_var_ub, decision_var,nparam)
! reading fun flag for optimization
read(24,*)runflag
if (runflag.eq.1) read(24, *) objective_function
close(24)

nres_level = my_reservoir(1)%nres_level
! total number of contraints 
ncons = nres + nuser + nres_level 
! storage reliability contraint for each reservoir
! user contract reliability contraint for each user
! restricition level probability contraints

Allocate(cons_global(ncons),value_net(nensem))

! The above line calculates the total number of constraints. It is better to represent equality constraint as inequality constraints
! It will ease the solver.
! At this point, rule curves are not incorporated as constraints. It could be modified if required.

Allocate(my_flow_set(nwatershed),parallel_track(nwatershed))

Do i = 1,nwatershed
	Allocate(my_flow_set(i)%controlled_flows(ntime), my_flow_set(i)%uncontrolled_flows(ntime,nensem))
end do

if(runflag == 1) then
	! Read parameters for FFSQP
	read(25,*)nf,mode,iprint,miter
	read(25,*)bigbnd,eps,epseqn,udelta 
	! model params to modify in the above file, if needed : 
	! nf : Number of objective fucntions 
	! mode : 110 [CBA - ref. ffsqp.f ] 
	!iprint : print level information 
	! miter : maximum Number of iteration 
	! bigbnd : plus infinity 
	! eps : stopping criterion that ensures a solution, the norm of the Newton direction vector is smaller than eps 
	! epseqn : tolerance of the violation of nonlinear equality constraints allowed by the user at an optimal solution 
	! udelta : perturbation size 
	close(25)
	! nineqn : Number of nonlinear inequality constraints
	nineqn = ncons 
	! nineq  : Number of inequality constraints     
	nineq = nineqn
	! neqn   : Number of nonlinear equality constraints 
	neqn = 0
	! neq    : Number of equality constraints 
	neq = 0
	! working space dimension allocation for fsqp 
	iwsize = 6*nparam + 8*max(1,nineq+neq)+7*max(1,nf)+30
	nwsize = 4*nparam*nparam + 5*max(1,nineq+neq)*nparam
	nwsize = nwsize + 3*max(1,nf)*nparam + 26*(nparam + max(1,nf))
	nwsize = nwsize + 45*max(1,nineq+neq) +100

	ncheck = nineq + neq
	! print *, 'ncheck = ',ncheck

	Allocate(iw(iwsize))
	Allocate(w(nwsize),g(ncheck),f(nf))

	! print *, 'iwsize = ',iwsize 
	! print *, 'nwsize = ',nwsize 
	print *, nparam

	ifinal = 0
	inform = 0
	index_cons = 1
	if (objective_function.eq."simple_benefits") then
		call simple_benefits(nparam, 1, decision_var, initial_bens)
		call FFSQP(nparam,nf,nineqn,nineq,neqn,neq,mode,iprint,miter,inform,bigbnd,eps,  &
				epsneq,udelta,decision_var_lb, decision_var_ub, decision_var,         &
				f,g,iw,iwsize,w,nwsize,simple_benefits,constr_res,grobfd,grcnfd)
		call simple_benefits(nparam, 1, decision_var, expect_bens)
		initial_bens = -initial_bens
		expect_bens = -expect_bens
	else if (objective_function.eq."max_release") then
		call max_release(nparam, 1, decision_var, initial_bens)
		call FFSQP(nparam,nf,nineqn,nineq,neqn,neq,mode,iprint,miter,inform,bigbnd,eps,  &
				epsneq,udelta,decision_var_lb, decision_var_ub, decision_var,         &
				f,g,iw,iwsize,w,nwsize,max_release,constr_res,grobfd,grcnfd)
		call max_release(nparam, 1, decision_var, expect_bens)
		initial_bens = -initial_bens
		expect_bens = -expect_bens
	else if (objective_function.eq."spill_objective") then
		call spill_objective(nparam, 1, decision_var, initial_bens)
		call FFSQP(nparam,nf,nineqn,nineq,neqn,neq,mode,iprint,miter,inform,bigbnd,eps,  &
				epsneq,udelta,decision_var_lb, decision_var_ub, decision_var,         &
				f,g,iw,iwsize,w,nwsize,spill_objective,constr_res,grobfd,grcnfd)
		call spill_objective(nparam, 1, decision_var, expect_bens)
	else
		call expected_benefits(nparam, 1, decision_var, initial_bens)
		call FFSQP(nparam,nf,nineqn,nineq,neqn,neq,mode,iprint,miter,inform,bigbnd,eps,  &
				epsneq,udelta,decision_var_lb, decision_var_ub, decision_var,         &
				f,g,iw,iwsize,w,nwsize,expected_benefits,constr_res,grobfd,grcnfd)
		call expected_benefits(nparam, 1, decision_var, expect_bens)
		initial_bens = -initial_bens
		expect_bens = -expect_bens
	end if
	write(38, "(F20.4)") expect_bens
end if 
ifinal = 1  
 
call constr_res(nparam,index_cons,decision_var,gcons)
write(*,"(A,1x,F12.2)") "Expected Benefits :", -expect_bens

write(33,50)(decision_var(k1),k1=1,icount_max)

close(31)
close(32)
close(33)
close(34)
close(35)
close(36)
close(37)
close(38)

21 FORMAT(I3,1x,A40)
22 FORMAT(I3,1x,A30) 
50 FORMAT(F10.3)

END 

! Following functions are for reading input files.

! subroutine to read watershed details.dat 
Subroutine read_watershed_details(my_watershed,nwatershed,ntime,nensem,nperiods)
USE Definitions
implicit doubleprecision(a-h,o-z)
character*40  file_name

TYPE(watershed) my_watershed(nwatershed)

read(11,*)iNum

read(11,20)my_watershed(iNum)%name

read(11,*)my_watershed(iNum)%ID, nchild,my_watershed(iNum)%drainage_area

my_watershed(iNum)%nchild = nchild

Allocate(my_watershed(iNum)%child_id(nchild))
Allocate(my_watershed(iNum)%child_type(nchild))


Do i = 1,nchild
	
	read(11,*)my_watershed(iNum)%child_type(i),my_watershed(iNum)%child_id(i)

END DO

Allocate(my_watershed(iNum)%natural_inflows(ntime*nperiods, nensem))

read(11,*)file_name

!print *, file_name
open(unit=40,file=trim(file_name))
Do k = 1,nensem
	Do i = i, nperiods
		read(40,*)(my_watershed(iNum)%natural_inflows(j,k), j=1+((i - 1)*ntime),ntime*i)
	End Do
End Do
close(40)

!Do k = 1,nensem

!	DO j = 1,ntime
      
!	  t1 = my_watershed(iNum)%natural_inflows(j,k)

!      if(j.eq.1.or.j.eq.2.or.j.eq.4.or.j.eq.6.or.j.eq.7.or.j.eq.9.or.j.eq.11)days = 31

!      if(j.eq.3.or.j.eq.5.or.j.eq.10.or.j.eq.12)days = 30

!      if(j.eq.8)days = 28
!		t1 = t1*days*24*3600/(10**6)

!		my_watershed(iNum)%natural_inflows(j,k) = t1

!	end do

!end do



20 FORMAT(A40)

RETURN

END
! subroutine to read natural flow details 
Subroutine read_nflow_details(my_natural_flow,nnatural_flow)
USE Definitions
implicit doubleprecision(a-h,o-z)

TYPE(Natural_flow) my_natural_flow(nnatural_flow)

read(12,*)iNum

read(12,20)my_natural_flow(iNum)%name

read(12,*)my_natural_flow(iNum)%start_type, my_natural_flow(iNum)%start_id
read(12,*)my_natural_flow(iNum)%end_type, my_natural_flow(iNum)%end_id


read(12,*)my_natural_flow(iNum)%minimum_discharge, my_natural_flow(iNum)%maximum_discharge

20 FORMAT(A40)

RETURN

END

! reads unit 13 : reservoir_details.dat 
Subroutine read_reservoir_details(my_reservoir,nres,ntime,nensem,nperiods)
USE Definitions
implicit doubleprecision(a-h,o-z)

TYPE(Reservoir) my_reservoir(nres)

read(13,*)iNum

read(13,20)my_reservoir(iNum)%name 
read(13,*)my_reservoir(iNum)%latitude, my_reservoir(iNum)%longitude
print *, my_reservoir(iNum)%name
read(13,*)my_reservoir(iNum)%elevation_max, my_reservoir(iNum)%elevation_min 
read(13,*)my_reservoir(iNum)%storage_max, my_reservoir(iNum)%storage_min, my_reservoir(iNum)%current_storage
read(13,*)my_reservoir(iNum)%elevation_storage_coeff(1),my_reservoir(iNum)%elevation_storage_coeff(2) &
		,my_reservoir(iNum)%elevation_storage_coeff(3)
read(13,*)my_reservoir(iNum)%storage_area_coeff(1),my_reservoir(iNum)%storage_area_coeff(2)

read(13,*)my_reservoir(iNum)%nspillway,my_reservoir(iNum)%Number_outlets
read(13,*)my_reservoir(iNum)%nres_level
read(13,*)my_reservoir(iNum)%nchild,my_reservoir(iNum)%nparent

n1 = my_reservoir(iNum)%nchild
n2 = my_reservoir(iNum)%nparent
n4 = my_reservoir(iNum)%nspillway

Allocate(my_reservoir(iNum)%child_id(n1),my_reservoir(iNum)%child_type(n1))
Allocate(my_reservoir(iNum)%parent_id(n2),my_reservoir(iNum)%parent_type(n2))
Allocate(my_reservoir(iNum)%spill_type(n4),my_reservoir(iNum)%crest_level(n4),my_reservoir(iNum)%discharge_max(n4))

Do i = 1,n4

	read(13,*)my_reservoir(iNum)%spill_type(i),my_reservoir(iNum)%crest_level(i), my_reservoir(iNum)%discharge_max(i)

End do

Do i = 1,my_reservoir(iNum)%nchild

	read(13,*)my_reservoir(iNum)%child_type(i), my_reservoir(iNum)%child_id(i)

end do

Do i = 1,my_reservoir(iNum)%nparent

	read(13,*)my_reservoir(iNum)%parent_type(i), my_reservoir(iNum)%parent_id(i)

end do

n3 = my_reservoir(iNum)%Number_outlets

Allocate(my_reservoir(iNum)%elevation_rvalve(n3),my_reservoir(iNum)%area_rvalve(n3))
Allocate(my_reservoir(iNum)%loss_coeff_max(n3),my_reservoir(iNum)%loss_coeff_min(n3))

Do i = 1, n3

	read(13,*)t1,t2,t3,t4,t5,t6
	my_reservoir(iNum)%elevation_rvalve(i) = t1
	my_reservoir(iNum)%area_rvalve(i)	   = t2
	my_reservoir(iNum)%loss_coeff_max(i)   = t3
	my_reservoir(iNum)%loss_coeff_max(i)   = t4
	my_reservoir(iNum)%target_storage      = t5
	my_reservoir(iNum)%storage_prob        = t6

End do

Allocate(my_reservoir(iNum)%rule_curve(ntime*nperiods),my_reservoir(iNum)%evaporation_rate(ntime*nperiods))

do j = 1, nperiods
	read(13,*)(my_reservoir(iNum)%rule_curve(i),i=1+((j-1)*ntime),ntime*j)
end do

do j = 1, nperiods
	read(13,*)(my_reservoir(iNum)%evaporation_rate(i),i=1+((j-1)*ntime),ntime*j)
end do

n4 = my_reservoir(iNum)%nres_level

Allocate(my_reservoir(iNum)%tar_restr_prob(n4))
	
read(13,*)(my_reservoir(iNum)%tar_restr_prob(i),i=1,n4)



20 FORMAT(A40)


RETURN

END

! reads user_details.dat 
Subroutine read_user_details(my_user,nuser,ntime,nensem,nperiods)
USE Definitions
implicit doubleprecision(a-h,o-z)

TYPE(User) my_user(nuser)

read(14,*)iNum

read(14,20)my_user(iNum)%name

read(14,*)my_user(iNum)%id, my_user(iNum)%user_type,my_user(iNum)%nchild, my_user(iNum)%nparent,my_user(iNum)%nres_level

n1 = my_user(iNum)%nchild
n2 = my_user(iNum)%nparent
nres_level = my_user(iNum)%nres_level


Allocate(my_user(iNum)%child_id(n1),my_user(iNum)%child_type(n1))
Allocate(my_user(iNum)%parent_id(n2),my_user(iNum)%parent_type(n2))
Allocate(my_user(iNum)%demand_fract(ntime*nperiods), my_user(iNum)%restr_fract(nres_level), my_user(iNum)%res_compensation(nres_level))

Do i = 1,my_user(iNum)%nchild

	read(14,*)my_user(iNum)%child_type(i), my_user(iNum)%child_id(i)

end do

Do i = 1,my_user(iNum)%nparent

	read(14,*)my_user(iNum)%parent_type(i), my_user(iNum)%parent_id(i)

end do

	read(14,*)t1,t2,t3,t4,t5,t6,t7

	my_user(iNum)%failure_prob      = (1-t1)
	my_user(iNum)%con_res_vol	   = t2
	my_user(iNum)%tariff           = t3
	my_user(iNum)%penalty		   = t4
	my_user(iNum)%minimum_release = t5
	my_user(iNum)%maximum_release = t6
	my_user(iNum)%penalty_compen  = t7
	do j = 1, nperiods
		read(14,*)(my_user(iNum)%demand_fract(i), i=1+((j-1)*ntime),ntime*j)
	end do
	read(14,*)(my_user(iNum)%restr_fract(i), i=1,nres_level)
	read(14,*)(my_user(iNum)%res_compensation(i), i=1,nres_level)
	
if(my_user(iNum)%user_type.eq.4)then

	read(14,*)t1,t2,t3,t4,t5,t6

	my_user(iNum)%max_discharge = t1
	my_user(iNum)%installed_capacity = t2
	my_user(iNum)%generator_efficiency = t3
	my_user(iNum)%storage_energy_coeff(1) = t4
	my_user(iNum)%storage_energy_coeff(2) = t5
	my_user(iNum)%unit_rate_energy = t6


        Allocate(my_user(iNum)%tail_elevation(ntime*nperiods))
!        Do k = 1,nensem
                
!        		read(14,*)(my_user(iNum)%tail_elevation(j1), j1=1,ntime)
!        END DO


end if

read(14,*)nlags

my_user(iNum)%nlags = nlags

Allocate(my_user(iNum)%ffraction(nlags))

read(14,*)(my_user(iNum)%ffraction(i),i=1,nlags)


20 FORMAT(A40)

RETURN

END

! reads node_details.dat 
Subroutine read_node_details(my_node,nfnode)
USE Definitions
implicit doubleprecision(a-h,o-z)

TYPE(flow_join_node) my_node(nfnode)

read(15,*)iNum

read(15,20)my_node(iNum)%name

read(15,*)my_node(iNum)%id, my_node(iNum)%nchild, my_node(iNum)%nparent

n1 = my_node(iNum)%nchild
n2 = my_node(iNum)%nparent


Allocate(my_node(iNum)%child_id(n1),my_node(iNum)%child_type(n1))
Allocate(my_node(iNum)%parent_id(n2),my_node(iNum)%parent_type(n2))

Do i = 1,my_node(iNum)%nchild

	read(15,*)my_node(iNum)%child_type(i), my_node(iNum)%child_id(i)

end do

Do i = 1,my_node(iNum)%nparent

	read(15,*)my_node(iNum)%parent_type(i), my_node(iNum)%parent_id(i)

end do

20 FORMAT(A40)

RETURN

END

Subroutine read_dir_inflows_details(my_dir_inflows,ndir_inflows)
USE Definitions
implicit doubleprecision(a-h,o-z)

TYPE(direct_inflows) my_dir_inflows(ndir_inflows)

read(16,*)iNum

read(16,20)my_dir_inflows(iNum)%name

read(16,*)my_dir_inflows(iNum)%start_type, my_dir_inflows(iNum)%start_id
read(16,*)my_dir_inflows(iNum)%end_type, my_dir_inflows(iNum)%end_id

read(16,*)my_dir_inflows(iNum)%minimum_discharge, my_dir_inflows(iNum)%maximum_discharge

20 FORMAT(A40)


RETURN

END

Subroutine read_ret_inflows_details(my_ret_inflows,nret_inflows)
USE Definitions
implicit doubleprecision(a-h,o-z)

TYPE(return_inflows) my_ret_inflows(nret_inflows)

read(17,*)iNum

read(17,20)my_ret_inflows(iNum)%name

read(17,*)my_ret_inflows(iNum)%start_type, my_ret_inflows(iNum)%start_id
read(17,*)my_ret_inflows(iNum)%end_type, my_ret_inflows(iNum)%end_id

read(17,*)my_ret_inflows(iNum)%minimum_discharge, my_ret_inflows(iNum)%maximum_discharge


20 FORMAT(A40)


RETURN

END

Subroutine read_diversions_details(my_diversions,ndiversion)
USE Definitions
implicit doubleprecision(a-h,o-z)

TYPE(diversion) my_diversions(ndiversion)

read(18,*)iNum

read(18,20)my_diversions(iNum)%name

read(18,*)my_diversions(iNum)%start_type, my_diversions(iNum)%start_id
read(18,*)my_diversions(iNum)%end_type, my_diversions(iNum)%end_id

read(18,*)my_diversions(iNum)%minimum_discharge, my_diversions(iNum)%maximum_discharge

20 FORMAT(A40)


RETURN

END

Subroutine read_spillflow_details(my_spill_flow,nspillflow)
USE Definitions
implicit doubleprecision(a-h,o-z)

TYPE(Spill_flow) my_spill_flow(nspillflow)

read(19,*)iNum

read(19,20)my_spill_flow(iNum)%name

read(19,*)my_spill_flow(iNum)%start_type, my_spill_flow(iNum)%start_id
read(19,*)my_spill_flow(iNum)%end_type, my_spill_flow(iNum)%end_id

read(19,*)my_spill_flow(iNum)%minimum_discharge, my_spill_flow(iNum)%maximum_discharge

20 FORMAT(A40)


RETURN

END

Subroutine read_interbasin_flow_details(my_interbasin_flow,ninterbasin_flow)
USE Definitions
implicit doubleprecision(a-h,o-z)

TYPE(Interbasin_flow) my_interbasin_flow(ninterbasin_flow)

read(20,*)iNum

read(20,20)my_interbasin_flow(iNum)%name

read(20,*)my_interbasin_flow(iNum)%start_type, my_interbasin_flow(iNum)%start_id
read(20,*)my_interbasin_flow(iNum)%end_type, my_interbasin_flow(iNum)%end_id

read(20,*)my_interbasin_flow(iNum)%minimum_discharge, my_interbasin_flow(iNum)%maximum_discharge

20 FORMAT(A40)


RETURN

END

Subroutine read_demandrelease_details(my_demand_release,ndemand_release)
USE Definitions
implicit doubleprecision(a-h,o-z)

TYPE(demand_release) my_demand_release(ndemand_release)

read(21,*)iNum

read(21,20)my_demand_release(iNum)%name

read(21,*)my_demand_release(iNum)%start_type, my_demand_release(iNum)%start_id
read(21,*)my_demand_release(iNum)%end_type, my_demand_release(iNum)%end_id

read(21,*)my_demand_release(iNum)%minimum_discharge, my_demand_release(iNum)%maximum_discharge

20 FORMAT(A40)


RETURN

END




Subroutine read_sink_details(my_sink,nsink)
USE Definitions
implicit doubleprecision(a-h,o-z)

TYPE(sink) my_sink(nsink)

read(22,*)iNum

read(22,*)my_sink(iNum)%name

read(22,*)my_sink(iNum)%id, my_sink(iNum)%nparent

n2 = my_sink(iNum)%nparent


Allocate(my_sink(iNum)%parent_id(n2),my_sink(iNum)%parent_type(n2))

Do i = 1,my_sink(iNum)%nparent

	read(22,*)my_sink(iNum)%parent_type(i), my_sink(iNum)%parent_id(i)

end do

read(22,*)my_sink(iNum)%max_storage

RETURN

END

Subroutine read_interbasin_details(my_interbasin,ninterbasin,ntime,nperiods)
USE Definitions
implicit doubleprecision(a-h,o-z)

TYPE(interbasin) my_interbasin(ninterbasin)

read(23,*)iNum

read(23,20)my_interbasin(iNum)%name
read(23,*) my_interbasin(iNum)%drainage_area

read(23,*)my_interbasin(iNum)%id, my_interbasin(iNum)%nchild

n1 = my_interbasin(iNum)%nchild


Allocate(my_interbasin(iNum)%child_id(n1),my_interbasin(iNum)%child_type(n1))

Do i = 1,my_interbasin(iNum)%nchild

	read(23,*)my_interbasin(iNum)%child_type(i), my_interbasin(iNum)%child_id(i)

end do

Allocate(my_interbasin(iNum)%average_flow(ntime*nperiods))
do i = 1, nperiods
	read(23,*)(my_interbasin(iNum)%average_flow(j), j=1+((i-1)*ntime),ntime*i)
end do

! outdated conversion
!Do j = 1,ntime

!	  t1 = my_interbasin(iNum)%average_flow(j)

!      if(j.eq.1.or.j.eq.2.or.j.eq.4.or.j.eq.6.or.j.eq.7.or.j.eq.9.or.j.eq.11)days = 31

!      if(j.eq.3.or.j.eq.5.or.j.eq.10.or.j.eq.12)days = 30

!      if(j.eq.8)days = 28

!		t1 = t1*days*24*3600/(10**6)

!		my_interbasin(iNum)%average_flow(j) = t1

!end do

20 FORMAT(A40)


RETURN

END
