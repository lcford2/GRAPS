!Program for multi-reservoir allocation
Program multireservoir
USE My_variables
implicit doubleprecision(a-h,o-z)
common /final_print/ifinal

external constr_res,expected_benefits,grobfd,grcnfd,ffsqp
!CHARACTER(*), PARAMETER :: inpath = "/input/"
!CHARACTER(*), PARAMETER :: outpath = "/output/"

double precision, allocatable	::	decision_var_lb(:),decision_var_ub(:), decision_var(:),w(:),g(:),f(:)
Integer,allocatable :: iw(:)

! subroutine declarations

open(unit=10, file = 'input.dat',				ACTION = 'READ', STATUS = 'OLD', ERR = 100)
open(unit=11, file = 'watershed_details.dat',	ACTION = 'READ', STATUS = 'OLD')
open(unit=12, file = 'nflow_details.dat',		ACTION = 'READ', STATUS = 'OLD')
open(unit=13, file = 'reservoir_details.dat',	ACTION = 'READ', STATUS = 'OLD')
open(unit=14, file = 'user_details.dat',		ACTION = 'READ', STATUS = 'OLD')
open(unit=15, file = 'node_details.dat',		ACTION = 'READ', STATUS = 'OLD')
open(unit=16, file = 'dir_flow_details.dat',	ACTION = 'READ', STATUS = 'OLD')
open(unit=17, file = 'ret_flow_details.dat',	ACTION = 'READ', STATUS = 'OLD')
open(unit=18, file = 'diversions_details.dat',	ACTION = 'READ', STATUS = 'OLD')
open(unit=19, file = 'spill_flow_details.dat',	ACTION = 'READ', STATUS = 'OLD')
open(unit=20, file = 'ibasin_flow_details.dat',	ACTION = 'READ', STATUS = 'OLD')
open(unit=21, file = 'demand_flow_details.dat',	ACTION = 'READ', STATUS = 'OLD')
open(unit=22, file = 'sink_details.dat',		ACTION = 'READ', STATUS = 'OLD')
open(unit=23, file = 'interbasin_details.dat',	ACTION = 'READ', STATUS = 'OLD')
open(unit=31, file = 'storage.out')
open(unit=54, file = 'hydro.out')
open(unit=24, file='release.out')
open(unit=28, file='flow_sets.out')
open(unit=29, file='node_flow.out')
open(unit=30, file='spill.out')
open(unit=44, file='deficit.out')
!  reading input.dat 
read(10,*)ntime,nensem
read(10,*)nwatershed,nnatural_flow,nres,nuser,nfnode,ndir_inflows,nret_inflows,ndiversion,nspill_flow,ninterbasin_flow,ndemand_release,nsink,ninterbasin

checksum = nwatershed+nnatural_flow+nres+nuser+nfnode+ndir_inflows+nret_inflows+ndiversion+nspill_flow+ninterbasin_flow+ndemand_release+nsink+ninterbasin

! ntime = number of time steps
! nres = number of reservoirs
! nuser = number of water users in the entire system
! nfnode = number of flow diversion/connection node
! nsink = number of sink points or netowrks ends
! nwatershed - NUmber of watershed originating natural flows
! ndir_inflows - NUmber of Direct Inflows released from the reservoir
! nret_flows = number of return flows (should not exceed nuser)
! ndiversion = NUmber of diversions towards environmental protection
! nspill_flow = Number of Spillways (should not exceed nres)
! nnatural_flow = Number of Natural flows from each watershed (should be equal to nwatershed)
! ninterbasin_flow = NUmber of Inter basin transfer points
! ndemand_release = NUmber of demand releases

! Allocate space for relevant data structures.

Allocate(my_user(nuser),my_reservoir(nres),my_node(nfnode),my_sink(nsink),my_watershed(nwatershed),my_interbasin(ninterbasin))

Allocate(my_dir_inflows(ndir_inflows),my_ret_inflows(nret_inflows),my_diversions(ndiversion))

Allocate(my_spill_flow(nspill_flow),my_natural_flow(nnatural_flow),my_interbasin_flow(ninterbasin_flow),&
    &my_demand_release(ndemand_release))

! Reads the system details from individual files until the connectivity given in input.dat comes to an end.
icount = 0
icount_max = checksum


DO WHILE (icount<icount_max)

! call each subroutines and read input files 

  read(10,21)itype,type_details
	if(itype.eq.1)call read_watershed_details(my_watershed,nwatershed,ntime,nensem)
	if(itype.eq.2)call read_nflow_details(my_natural_flow,nnatural_flow)
	if(itype.eq.3)call read_reservoir_details(my_reservoir,nres,ntime,nensem)
	if(itype.eq.4)call read_user_details(my_user,nuser,ntime, nensem)
  if(itype.eq.5)call read_node_details(my_node,nfnode)
	if(itype.eq.6)call read_dir_inflows_details(my_dir_inflows,ndir_inflows)
	if(itype.eq.7)call read_ret_inflows_details(my_ret_inflows,nret_inflows)
	if(itype.eq.8)call read_diversions_details(my_diversions,ndiversion)
	if(itype.eq.9)call read_spillflow_details(my_spill_flow,nspill_flow)
	if(itype.eq.10)call read_interbasin_flow_details(my_interbasin_flow,ninterbasin_flow)
	if(itype.eq.11)call read_demandrelease_details(my_demand_release,ndemand_release)
	if(itype.eq.12)call read_sink_details(my_sink,nsink)
	if(itype.eq.13)call read_interbasin_details(my_interbasin,ninterbasin,ntime)

icount = icount +1 

END DO

!print *, 'wwww'
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

nparam = (nuser*ntime) ! each reservoir should be given a minimum of one user for downstream release and it should be defined in the user_details.dat
nsimul_block = nres + nfnode
ntotal_vertices = nwatershed + nres + nuser + nfnode + ninterbasin_flow + nsink 

Allocate(decision_var_lb(nparam), decision_var_ub(nparam), decision_var(nparam))

Allocate(temp_decision_var(nparam),isimul_status(nsimul_block),my_network_order(nsimul_block))

Allocate(searched_vertices(ntotal_vertices))

! Initialize temporary decision variable

call array_ini(nparam,temp_decision_var,-10.0d0)

call solution_path()

print *, (k1, my_network_order(k1)%order_type, my_network_order(k1)%order_id, k1=1, nsimul_block)

! Extract lower and upper bounds of decision variables.

call constr_decision_extract(decision_var_lb,decision_var_ub, decision_var,nparam)
open(unit=60,file='runflag.dat',ACTION = 'READ', STATUS = 'OLD')
read(60,*)runflag
close(60)

if(runflag == 1)then 
! Read parameters for FFSQP

open(unit =40, file='model_para.dat',ACTION = 'READ', STATUS = 'OLD')

read(40,*)nf,mode,iprint,miter
read(40,*)bigbnd,eps,epseqn,udelta 
! model params to modify in the above file, if needed : 
! nf : number of objective fucntions 
! mode : 110 [CBA - ref. ffsqp.f ] 
!iprint : print level information 
! miter : maximum number of iteration 
! bigbnd : plus infinity 
! eps : stopping criterion that ensures a solution, the norm of the Newton direction vector is smaller than eps 
! epseqn : tolerance of the violation of nonlinear equality constraints allowed by the user at an optimal solution 
! udelta : perturbation size 

close(40)
end if 

nres_level = my_reservoir(1)%nres_level
! total number of constraints 
ncons = nuser + nres_level + nres  

Allocate(cons_global(ncons),value_net(nensem))

! The above line calculates the total number of constraints. It is better to represent equality constraint as inequality constraints
! It will ease the solver.
! nuser - represents reliability constraint for each user
! nres_level - represents the total number of restricion level constraints
! nres - represents the end of the year storage constraints
! AT this point, rule curves are not incorporated as constraints. It could be modified if required.

Allocate(my_flow_set(nwatershed),parallel_track(nwatershed))


Do i = 1,nwatershed

		Allocate(my_flow_set(i)%controlled_flows(ntime), my_flow_set(i)%uncontrolled_flows(ntime,nensem))

end do
if(runflag == 1) then
! nineqn : number of nonlinear inequality constraints
nineqn = ncons 
! nineq  : number of inequality constraints     
nineq = nineqn
! neqn   : number of nonlinear equality constraints 
neqn = 0
! neq    : number of equality constraints 
neq = 0
! working space dimension allocation for fsqp 
iwsize = 6*nparam + 8*max(1,nineq+neq)+7*max(1,nf)+30
nwsize = 4*nparam*nparam + 5*max(1,nineq+neq)*nparam
nwsize = nwsize + 3*max(1,nf)*nparam + 26*(nparam + max(1,nf))
nwsize = nwsize + 45*max(1,nineq+neq) +100

ncheck = nineq + neq
print *, 'ncheck = ',ncheck

Allocate(iw(iwsize))
Allocate(w(nwsize),g(ncheck),f(nf))

print *, 'iwsize = ',iwsize 
print *, 'nwsize = ',nwsize 
print *, nparam

ifinal = 0
inform = 0
index_cons = 1
print *, "entering FSQP..........."
call FFSQP(nparam,nf,nineqn,nineq,neqn,neq,mode,iprint,miter,inform,bigbnd,eps,  &
           epsneq,udelta,decision_var_lb, decision_var_ub, decision_var,         &
          f,g,iw,iwsize,w,nwsize,expected_benefits,constr_res,grobfd,grcnfd)

print *, "Completed running FSQP............."
end if 
ifinal = 1  
 
call constr_res(nparam,index_cons,decision_var,gcons)
print *, "ifinal = ",ifinal
! end running optimization or simulation
print *, "exiting..........."
!write(24,*)'Constraints :'
!write(24,*)(cons_global(k1),k1=1,ncons)
!write(24,*)'Objective function'
!write(24,*)(value_net(k1),k1=1,nensem)
!write(24,*) 'Optimised decision variables'
!icount = 1
icount_max = nuser*ntime
!do while (icount <= icount_max)
write(24,50)(decision_var(k1),k1=1,icount_max)
!icount = icount + 1
!end do 

close(24)
21  FORMAT(I3,1x,A40)
22 FORMAT(I3,1x,A30) 
50 FORMAT(F10.3)
100	STOP 
! end of multireservoir 
 	END 

! Following functions are for reading input files.

! subroutine to read watershed details.dat 
Subroutine read_watershed_details(my_watershed,nwatershed,ntime, nensem)
USE Definitions
implicit doubleprecision(a-h,o-z)
character*40  file_name

TYPE(watershed) my_watershed(nwatershed)

read(11,*)inum

read(11,20)my_watershed(inum)%name

read(11,*)my_watershed(inum)%ID, nchild,my_watershed(inum)%drainage_area

my_watershed(inum)%nchild = nchild

Allocate(my_watershed(inum)%child_id(nchild))
Allocate(my_watershed(inum)%child_type(nchild))


Do i = 1,nchild
	
	read(11,*)my_watershed(inum)%child_type(i),my_watershed(inum)%child_id(i)

END DO

Allocate(my_watershed(inum)%natural_inflows(ntime, nensem))

read(11,*)file_name

!print *, file_name
open(unit=40,file=trim(file_name))

Do k = 1,nensem
                
		read(40,*)(my_watershed(inum)%natural_inflows(j,k), j=1,ntime)


END DO

close(40)

!Do k = 1,nensem

!	DO j = 1,ntime
      
!	  t1 = my_watershed(inum)%natural_inflows(j,k)

!      if(j.eq.1.or.j.eq.2.or.j.eq.4.or.j.eq.6.or.j.eq.7.or.j.eq.9.or.j.eq.11)days = 31

!      if(j.eq.3.or.j.eq.5.or.j.eq.10.or.j.eq.12)days = 30

!      if(j.eq.8)days = 28
!		t1 = t1*days*24*3600/(10**6)

!		my_watershed(inum)%natural_inflows(j,k) = t1

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

read(12,*)inum

read(12,20)my_natural_flow(inum)%name

read(12,*)my_natural_flow(inum)%start_type, my_natural_flow(inum)%start_id
read(12,*)my_natural_flow(inum)%end_type, my_natural_flow(inum)%end_id


read(12,*)my_natural_flow(inum)%minimum_discharge, my_natural_flow(inum)%maximum_discharge

20 FORMAT(A40)

RETURN

END

! reads unit 13 : reservoir_details.dat 
Subroutine read_reservoir_details(my_reservoir,nres,ntime,nensem)
USE Definitions
implicit doubleprecision(a-h,o-z)

TYPE(Reservoir) my_reservoir(nres)

read(13,*)inum

read(13,20)my_reservoir(inum)%name 
read(13,*)my_reservoir(inum)%latitude, my_reservoir(inum)%longitude
print *, my_reservoir(inum)%name
read(13,*)my_reservoir(inum)%elevation_max, my_reservoir(inum)%elevation_min 
read(13,*)my_reservoir(inum)%storage_max, my_reservoir(inum)%storage_min, my_reservoir(inum)%current_storage
read(13,*)my_reservoir(inum)%elevation_storage_coeff(1),my_reservoir(inum)%elevation_storage_coeff(2) &
		,my_reservoir(inum)%elevation_storage_coeff(3)
read(13,*)my_reservoir(inum)%storage_area_coeff(1),my_reservoir(inum)%storage_area_coeff(2)

read(13,*)my_reservoir(inum)%nspillway,my_reservoir(inum)%number_outlets
read(13,*)my_reservoir(inum)%nres_level
read(13,*)my_reservoir(inum)%nchild,my_reservoir(inum)%nparent

n1 = my_reservoir(inum)%nchild
n2 = my_reservoir(inum)%nparent
n4 = my_reservoir(inum)%nspillway

Allocate(my_reservoir(inum)%child_id(n1),my_reservoir(inum)%child_type(n1))
Allocate(my_reservoir(inum)%parent_id(n2),my_reservoir(inum)%parent_type(n2))
Allocate(my_reservoir(inum)%spill_type(n4),my_reservoir(inum)%crest_level(n4),my_reservoir(inum)%discharge_max(n4))

Do i = 1,n4

	read(13,*)my_reservoir(inum)%spill_type(i),my_reservoir(inum)%crest_level(i), my_reservoir(inum)%discharge_max(i)

End do

Do i = 1,my_reservoir(inum)%nchild

	read(13,*)my_reservoir(inum)%child_type(i), my_reservoir(inum)%child_id(i)

end do

Do i = 1,my_reservoir(inum)%nparent

	read(13,*)my_reservoir(inum)%parent_type(i), my_reservoir(inum)%parent_id(i)

end do

n3 = my_reservoir(inum)%number_outlets

Allocate(my_reservoir(inum)%elevation_rvalve(n3),my_reservoir(inum)%area_rvalve(n3))
Allocate(my_reservoir(inum)%loss_coeff_max(n3),my_reservoir(inum)%loss_coeff_min(n3))

Do i = 1, n3

	read(13,*)t1,t2,t3,t4,t5,t6
	my_reservoir(inum)%elevation_rvalve(i) = t1
	my_reservoir(inum)%area_rvalve(i)	   = t2
	my_reservoir(inum)%loss_coeff_max(i)   = t3
	my_reservoir(inum)%loss_coeff_max(i)   = t4
	my_reservoir(inum)%target_storage      = t5
	my_reservoir(inum)%storage_prob        = t6

End do

Allocate(my_reservoir(inum)%rule_curve(ntime),my_reservoir(inum)%evaporation_rate(ntime))

read(13,*)(my_reservoir(inum)%rule_curve(i),i=1,ntime)
read(13,*)(my_reservoir(inum)%evaporation_rate(i),i=1,ntime)

n4 = my_reservoir(inum)%nres_level

Allocate(my_reservoir(inum)%tar_restr_prob(n4))
	
read(13,*)(my_reservoir(inum)%tar_restr_prob(i),i=1,n4)



20 FORMAT(A40)


RETURN

END

! reads user_details.dat 
Subroutine read_user_details(my_user,nuser,ntime, nensem)
USE Definitions
implicit doubleprecision(a-h,o-z)

TYPE(User) my_user(nuser)

read(14,*)inum

read(14,20)my_user(inum)%name

read(14,*)my_user(inum)%id, my_user(inum)%user_type,my_user(inum)%nchild, my_user(inum)%nparent,my_user(inum)%nres_level

n1 = my_user(inum)%nchild
n2 = my_user(inum)%nparent
nres_level = my_user(inum)%nres_level


Allocate(my_user(inum)%child_id(n1),my_user(inum)%child_type(n1))
Allocate(my_user(inum)%parent_id(n2),my_user(inum)%parent_type(n2))
Allocate(my_user(inum)%demand_fract(ntime), my_user(inum)%restr_fract(nres_level), my_user(inum)%res_compensation(nres_level))

Do i = 1,my_user(inum)%nchild

	read(14,*)my_user(inum)%child_type(i), my_user(inum)%child_id(i)

end do

Do i = 1,my_user(inum)%nparent

	read(14,*)my_user(inum)%parent_type(i), my_user(inum)%parent_id(i)

end do

	read(14,*)t1,t2,t3,t4,t5,t6,t7

	my_user(inum)%failure_prob      = (1-t1)
	my_user(inum)%con_res_vol	   = t2
	my_user(inum)%tariff           = t3
	my_user(inum)%penalty		   = t4
	my_user(inum)%minimum_release = t5
	my_user(inum)%maximum_release = t6
	my_user(inum)%penalty_compen  = t7
	read(14,*)(my_user(inum)%demand_fract(i), i=1,ntime)
	read(14,*)(my_user(inum)%restr_fract(i), i=1,nres_level)
	read(14,*)(my_user(inum)%res_compensation(i), i=1,nres_level)
	
if(my_user(inum)%user_type.eq.4)then

	read(14,*)t1,t2,t3,t4,t5,t6

	my_user(inum)%max_discharge = t1
	my_user(inum)%installed_capacity = t2
	my_user(inum)%generator_efficiency = t3
	my_user(inum)%storage_energy_coeff(1) = t4
	my_user(inum)%storage_energy_coeff(2) = t5
	my_user(inum)%unit_rate_energy = t6


        Allocate(my_user(inum)%tail_elevation(ntime))
!        Do k = 1,nensem
                
!        		read(14,*)(my_user(inum)%tail_elevation(j1), j1=1,ntime)
!        END DO


end if

read(14,*)nlags

my_user(inum)%nlags = nlags

Allocate(my_user(inum)%ffraction(nlags))

read(14,*)(my_user(inum)%ffraction(i),i=1,nlags)


20 FORMAT(A40)

RETURN

END

! reads node_details.dat 
Subroutine read_node_details(my_node,nfnode)
USE Definitions
implicit doubleprecision(a-h,o-z)

TYPE(flow_join_node) my_node(nfnode)

read(15,*)inum

read(15,20)my_node(inum)%name

read(15,*)my_node(inum)%id, my_node(inum)%nchild, my_node(inum)%nparent

n1 = my_node(inum)%nchild
n2 = my_node(inum)%nparent


Allocate(my_node(inum)%child_id(n1),my_node(inum)%child_type(n1))
Allocate(my_node(inum)%parent_id(n2),my_node(inum)%parent_type(n2))

Do i = 1,my_node(inum)%nchild

	read(15,*)my_node(inum)%child_type(i), my_node(inum)%child_id(i)

end do

Do i = 1,my_node(inum)%nparent

	read(15,*)my_node(inum)%parent_type(i), my_node(inum)%parent_id(i)

end do

20 FORMAT(A40)

RETURN

END

Subroutine read_dir_inflows_details(my_dir_inflows,ndir_inflows)
USE Definitions
implicit doubleprecision(a-h,o-z)

TYPE(direct_inflows) my_dir_inflows(ndir_inflows)

read(16,*)inum

read(16,20)my_dir_inflows(inum)%name

read(16,*)my_dir_inflows(inum)%start_type, my_dir_inflows(inum)%start_id
read(16,*)my_dir_inflows(inum)%end_type, my_dir_inflows(inum)%end_id

read(16,*)my_dir_inflows(inum)%minimum_discharge, my_dir_inflows(inum)%maximum_discharge

20 FORMAT(A40)


RETURN

END

Subroutine read_ret_inflows_details(my_ret_inflows,nret_inflows)
USE Definitions
implicit doubleprecision(a-h,o-z)

TYPE(return_inflows) my_ret_inflows(nret_inflows)

read(17,*)inum

read(17,20)my_ret_inflows(inum)%name

read(17,*)my_ret_inflows(inum)%start_type, my_ret_inflows(inum)%start_id
read(17,*)my_ret_inflows(inum)%end_type, my_ret_inflows(inum)%end_id

read(17,*)my_ret_inflows(inum)%minimum_discharge, my_ret_inflows(inum)%maximum_discharge


20 FORMAT(A40)


RETURN

END

Subroutine read_diversions_details(my_diversions,ndiversion)
USE Definitions
implicit doubleprecision(a-h,o-z)

TYPE(diversion) my_diversions(ndiversion)

read(18,*)inum

read(18,20)my_diversions(inum)%name

read(18,*)my_diversions(inum)%start_type, my_diversions(inum)%start_id
read(18,*)my_diversions(inum)%end_type, my_diversions(inum)%end_id

read(18,*)my_diversions(inum)%minimum_discharge, my_diversions(inum)%maximum_discharge

20 FORMAT(A40)


RETURN

END

Subroutine read_spillflow_details(my_spill_flow,nspillflow)
USE Definitions
implicit doubleprecision(a-h,o-z)

TYPE(Spill_flow) my_spill_flow(nspillflow)

read(19,*)inum

read(19,20)my_spill_flow(inum)%name

read(19,*)my_spill_flow(inum)%start_type, my_spill_flow(inum)%start_id
read(19,*)my_spill_flow(inum)%end_type, my_spill_flow(inum)%end_id

read(19,*)my_spill_flow(inum)%minimum_discharge, my_spill_flow(inum)%maximum_discharge

20 FORMAT(A40)


RETURN

END

Subroutine read_interbasin_flow_details(my_interbasin_flow,ninterbasin_flow)
USE Definitions
implicit doubleprecision(a-h,o-z)

TYPE(Interbasin_flow) my_interbasin_flow(ninterbasin_flow)

read(20,*)inum

read(20,20)my_interbasin_flow(inum)%name

read(20,*)my_interbasin_flow(inum)%start_type, my_interbasin_flow(inum)%start_id
read(20,*)my_interbasin_flow(inum)%end_type, my_interbasin_flow(inum)%end_id

read(20,*)my_interbasin_flow(inum)%minimum_discharge, my_interbasin_flow(inum)%maximum_discharge

20 FORMAT(A40)


RETURN

END

Subroutine read_demandrelease_details(my_demand_release,ndemand_release)
USE Definitions
implicit doubleprecision(a-h,o-z)

TYPE(demand_release) my_demand_release(ndemand_release)

read(21,*)inum

read(21,20)my_demand_release(inum)%name

read(21,*)my_demand_release(inum)%start_type, my_demand_release(inum)%start_id
read(21,*)my_demand_release(inum)%end_type, my_demand_release(inum)%end_id

read(21,*)my_demand_release(inum)%minimum_discharge, my_demand_release(inum)%maximum_discharge

20 FORMAT(A40)


RETURN

END




Subroutine read_sink_details(my_sink,nsink)
USE Definitions
implicit doubleprecision(a-h,o-z)

TYPE(sink) my_sink(nsink)

read(22,*)inum

read(22,*)my_sink(inum)%name

read(22,*)my_sink(inum)%id, my_sink(inum)%nparent

n2 = my_sink(inum)%nparent


Allocate(my_sink(inum)%parent_id(n2),my_sink(inum)%parent_type(n2))

Do i = 1,my_sink(inum)%nparent

	read(22,*)my_sink(inum)%parent_type(i), my_sink(inum)%parent_id(i)

end do

read(22,*)my_sink(inum)%max_storage

RETURN

END

Subroutine read_interbasin_details(my_interbasin,ninterbasin,ntime)
USE Definitions
implicit doubleprecision(a-h,o-z)

TYPE(interbasin) my_interbasin(ninterbasin)

read(23,*)inum

read(23,20)my_interbasin(inum)%name
read(23,*) my_interbasin(inum)%drainage_area

read(23,*)my_interbasin(inum)%id, my_interbasin(inum)%nchild

n1 = my_interbasin(inum)%nchild


Allocate(my_interbasin(inum)%child_id(n1),my_interbasin(inum)%child_type(n1))

Do i = 1,my_interbasin(inum)%nchild

	read(23,*)my_interbasin(inum)%child_type(i), my_interbasin(inum)%child_id(i)

end do

Allocate(my_interbasin(inum)%average_flow(ntime))

read(23,*)(my_interbasin(inum)%average_flow(j), j=1,ntime)

! outdated conversion
!Do j = 1,ntime

!	  t1 = my_interbasin(inum)%average_flow(j)

!      if(j.eq.1.or.j.eq.2.or.j.eq.4.or.j.eq.6.or.j.eq.7.or.j.eq.9.or.j.eq.11)days = 31

!      if(j.eq.3.or.j.eq.5.or.j.eq.10.or.j.eq.12)days = 30

!      if(j.eq.8)days = 28

!		t1 = t1*days*24*3600/(10**6)

!		my_interbasin(inum)%average_flow(j) = t1

!end do

20 FORMAT(A40)


RETURN

END
