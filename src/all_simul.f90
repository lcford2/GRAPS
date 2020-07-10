subroutine constr_res(nparam,index_cons,decision_var,gcons)
use My_variables
implicit doubleprecision(a-h,o-z)
common /final_print/ifinal

Real*8 decision_var(nparam),total_deficit(nensem),gcons
! integer ifinal
! check whether the current value of iteration is same as that of previous iteration
! if so , just return the current value of the constraint using the variable index_cons.


ivar_status = icheck_var_status(nparam,decision_var,temp_decision_var)
if (ifinal.eq.1) ivar_status =1

if(ivar_status.eq.0)then
	gcons = cons_global(index_cons)
    return
end if

! if the decision_variable value has changed then perform simulation
nsimul_block = nres + nfnode


itemp_type = 0
itemp_id = 0

iflow_set = 0
icurrent_id = 0
icurrent_type = 0


call array_ini(nensem,total_deficit,0.0d0)
call array_ini(nensem,value_net,0.0d0)
call array_ini(ncons,cons_global,0.0d0)

! Initialize controlled and uncontrolled flows for each watershed

do i = 1,nwatershed
	call array_ini(ntime,my_flow_set(i)%controlled_flows,0.0d0)
	do j = 1,nensem
		call array_ini(ntime,my_flow_set(i)%uncontrolled_flows(1,j),0.0d0)
	end do
	parallel_track(i)%order_type = 0
	parallel_track(i)%order_id = 0
end do

icount = 0

! Loop for simulation at each junction node and reservoir
do i = 1, nsimul_block
	iprev_id = icurrent_id
	iprev_type = icurrent_type

	icurrent_type = my_network_order(i)%order_type
	icurrent_id = my_network_order(i)%order_id

	if(icurrent_type.eq.5)then	
	    call node_simul_module(icurrent_type,icurrent_id,iprev_id,iprev_type,nparam,decision_var)
	else if(icurrent_type.eq.3) then
		call reservoir_simul_module(icurrent_type,icurrent_id, iprev_id, &
	   			iprev_type,nparam,decision_var,total_deficit,nend_cons)
		icount = icount +1
		ensem = nensem
		icount = icount +1
		ensem = nensem
		! Calculation for end of time steps storage constraints
		! If it is only a single ensemble member, don't deal with probabilities.
		if (nensem.eq.1) then
			cons_global(icount) = my_reservoir(icurrent_id)%target_storage - my_reservoir(icurrent_id)%eot_storage
		else
			cons_global(icount) = (real(nend_cons)/real(ensem)) - my_reservoir(icurrent_id)%storage_prob
		end if
	end if	
end do

call deficit_splitter(total_deficit,decision_var,nparam)

if(ifinal.eq.1) print *, index_cons
gcons = cons_global(index_cons)
return
end

! subroutine for reservoir simulation
subroutine reservoir_simul_module(icurrent_type,icurrent_id,iprev_id, &
	   iprev_type, nparam,decision_var,total_deficit,nend_cons)
use My_variables
implicit doubleprecision(a-h,o-z)
common /final_print/ifinal

double precision decision_var(nparam),q(ntime),simul_deficit(ntime),simul_evapo(ntime)
double precision simul_stor(ntime), simul_spill(ntime),rate_area(ntime),release(ntime),act_release(ntime)
double precision total_deficit(nensem)

nparent = my_reservoir(icurrent_id)%nparent

ijump = 0
iadd = 0

do i1 = 1,nparent
	if((iprev_type.eq.my_reservoir(icurrent_id)%parent_type(i1)).and. &
		(iprev_id.eq.my_reservoir(icurrent_id)%parent_id(i1)))ijump = 1
end do

if(ijump.eq.0)then	
	if((iprev_type.ne.0).and.(iprev_id.ne.0))then
		parallel_track(iflow_set)%order_id = iprev_id
		parallel_track(iflow_set)%order_type = iprev_type
	end if
	iflow_set = iflow_set +1        
end if 

do i = 1,nensem
	do j = 1,ntime
		my_flow_set(iflow_set)%uncontrolled_flows(j,i) = 0
	end do
end do

inodeparent = 0

do i1 = 1,nparent
	iparent_type = my_reservoir(icurrent_id)%parent_type(i1)
	if (iparent_type == 5)  then 
		inodeparent = 1
	end if
end do


if (inodeparent.ne.1) then
	do j = 1,ntime
		my_flow_set(iflow_set)%controlled_flows(j) = 0
	end do
end if

! Loop for adding controlled and uncontrolled flow from parents
do i1 = 1,nparent
	iparent_type = my_reservoir(icurrent_id)%parent_type(i1)
	iparent_id   = my_reservoir(icurrent_id)%parent_id(i1)
	! Prepare the inflow sets

	iadd = 0

	! If parent is a watershed
	if ((iparent_type.eq.1).or.(iparent_type.eq.3)) then
		call add_uncontrolled_flows(iparent_type,iparent_id,decision_var,nparam)
	end if

	! If parent is another reservoir
	if(iparent_type.eq.3) then
		call add_controlled_flows(iparent_type,iparent_id, decision_var,nparam)
		! call add_spill_flows(iparent_id)
	end if
	
	! If parent is a junction node
	if(iparent_type.eq.5) then
		do j1 = 1,nwatershed
			if((parallel_track(j1)%order_type.eq.iparent_type).and.(parallel_track(j1)%order_id.eq.iparent_id)) then 
				call add_all_flows(j1,iparent_type,iparent_id,decision_var,nparam)
			end if
		end do
		if((iparent_type.ne.iprev_type).and.(iparent_id.ne.iprev_id)) then
			call add_controlled_flows(iparent_type,iparent_id, decision_var,nparam)
		end if
	end if
	! if parent is a user or a interbasin transfer
	if((iparent_type.eq.4).or.(iparent_type.eq.13)) then
		call add_controlled_flows(iparent_type,iparent_id, decision_var,nparam)
	end if
end do

! Prepare the outflow sets
call array_ini(ntime,release, 0.0d0)
call array_ini(ntime,act_release, 0.0d0)

nchild = my_reservoir(icurrent_id)%nchild
       
! Loop for adding the releases to child nodes
do i1 = 1,nchild
	ichild_type = my_reservoir(icurrent_id)%child_type(i1)
	ichild_id   = my_reservoir(icurrent_id)%child_id(i1)
	call calculate_outflows(ichild_type,ichild_id,icurrent_id, icurrent_type, release,decision_var,nparam)
end do

storage_max = my_reservoir(icurrent_id)%storage_max
storage_ini = my_reservoir(icurrent_id)%current_storage

do j = 1,ntime
	rate_area(j) = my_reservoir(icurrent_id)%evaporation_rate(j)
end do


! Variable for counting the number of time not meeting target storage constraints
nend_cons = 0

! Loop for reservoir simulation
do i = 1,nensem
	do j = 1,ntime
		q(j) = my_flow_set(iflow_set)%uncontrolled_flows(j,i) + &
		       my_flow_set (iflow_set)%controlled_flows(j)
	end do
	call reser_simul(icurrent_id, q, ntime,release,storage_max, storage_ini,rate_area, &
					 simul_stor,simul_spill, simul_evapo, simul_deficit,iflag)
	do j=1,ntime
		spill_values((icurrent_id - 1) * ntime + j) = simul_spill(j)
		my_reservoir(icurrent_id)%spill(j, i) = simul_spill(j)
		act_release(j) = release(j)+simul_spill(j)-simul_deficit(j)
		if(act_release(j).le.0.0)act_release(j) = 0.0
	end do
	my_reservoir(icurrent_id)%eot_storage = simul_stor(ntime)

  	if (ifinal.eq.1) then
    	write(31,15) my_reservoir(icurrent_id)%name, (simul_stor(j), j=1,ntime)        
			15 format(a,2x, <ntime>(2x,F15.3))
		write(34,15) my_reservoir(icurrent_id)%name, (act_release(j), j=1,ntime)
		write(37,15) my_reservoir(icurrent_id)%name, (simul_deficit(j), j=1,ntime)
		write(36,15) my_reservoir(icurrent_id)%name, (simul_spill(j), j=1, ntime)
		if (nensem.eq.1)  then
			print *, my_reservoir(icurrent_id)%name, 'Final Storage:',simul_stor(ntime)
			print *, my_reservoir(icurrent_id)%name, 'Total Spill:', sum(simul_spill)
		end if
	end if  
	
	do i1 = 1,nchild
		ichild_type = my_reservoir(icurrent_id)%child_type(i1)
		ichild_id   = my_reservoir(icurrent_id)%child_id(i1)
		if (ichild_type.eq.4) then
			if (my_user(ichild_id)%user_type.eq.4) then
				call hydropower(ichild_type,ichild_id,icurrent_id, icurrent_type, simul_stor, release, nensem)
			end if 
		end if
	end do

	if(iflag.eq.1)then
		do k=1,ntime
			total_deficit(i) = total_deficit(i)+simul_deficit(k)
		end do
	end if

	if(simul_stor(ntime).lt.my_reservoir(icurrent_id)%target_storage) then
		nend_cons = nend_cons + 1
	end if
end do

if (ifinal.eq.1) then
	if(nensem.gt.1.0)	print *, my_reservoir(icurrent_id)%name, "Reliability =", 100*(1-(real(nend_cons)/real(nensem))),'%'
end if

return
end


! subroutine for hydropwer calculation
subroutine hydropower(iblock_type,iblock_id,icurrent_id, icurrent_type,simul_stor,release,nensemble)
use My_variables
implicit doubleprecision(a-h,o-z)
common /final_print/ifinal
double precision simul_hydropower(ntime), simul_stor(ntime), release(ntime)
double precision head_elevation(ntime)
double precision Unit_Conv
		! output in MWh/month
        !Unit_Conv = 62.4*43560*0.746/3600/550 ! Amir's Unit Conv =  62.4*43560*1000/0.74/3600/1000/1000  ! Yi's Unit_Conv = 62.5*1.356*1000/2/31/1000000
        Unit_Conv = 9.81*1000/3600 ! [9.81 KN/m^3]*[10^6m^3/hm^3]/[3600 KJ/KWh]*[1000 KWh/MWh] flow in hm3 per month, head in m -Lucas
        stor_uni = my_reservoir(icurrent_id)%current_storage
        do j = 1, ntime
            !Calculating head Elevation
                
            tmp = (stor_uni + simul_stor(j)) /2
            head_elevation(j)=my_reservoir(icurrent_id)%elevation_storage_coeff(1)*(tmp**2)+&
                        my_reservoir(icurrent_id)%elevation_storage_coeff(2)*tmp+ &
                        my_reservoir(icurrent_id)%elevation_storage_coeff(3)

            !Calculating Simulated Hydropower
            simul_hydropower(j) =  my_user(iblock_id)%generator_efficiency*release(j)*& 
            ( head_elevation(j) - my_user(iblock_id)%tail_elevation(j))*Unit_Conv

                

            if  (simul_hydropower(j) .LT. 0.0)then 
                    simul_hydropower(j) = 0.0
            end if 
     
            stor_uni = simul_stor(j)
        end do

   if (ifinal.eq.1) then 
       write(32,19) my_user(icurrent_id)%name, (simul_hydropower(j), j=1,ntime) 
        19 format(a,2x, <ntime>(2x,F16.3))
        end if 
return
end

! subroutine for simulation at each junction node

subroutine node_simul_module(icurrent_type,icurrent_id,iprev_id, &
	   iprev_type, nparam,decision_var)
use My_variables
implicit doubleprecision(a-h,o-z)
common /final_print/ifinal
double precision decision_var(nparam),release(ntime),inflow(ntime)
double precision user_release(ntime), remaining_flow(ntime, nensem), loss_factor
character*40 output_string
integer cons_loc

nparent = my_node(icurrent_id)%nparent

ijump = 0
iadd = 0

do i1 = 1,nparent		
	if((iprev_type.eq.my_node(icurrent_id)%parent_type(i1)).and. &
		(iprev_id.eq.my_node(icurrent_id)%parent_id(i1)))ijump = 1
end do

if(ijump.eq.0)then
	if((iprev_type.ne.0).and.(iprev_id.ne.0))then
		parallel_track(iflow_set)%order_id = iprev_id
		parallel_track(iflow_set)%order_type = iprev_type
	end if
	iflow_set = iflow_set + 1
end if 

if (iflow_set.ne.1) then
	do j = 1,ntime
		my_flow_set (iflow_set-1)%controlled_flows(j) = 0.0
	end do
end if 

do j = 1,ntime
	my_flow_set (iflow_set)%controlled_flows(j) = 0.0
end do

! Add controlled and uncontrolled flows from parents
do i1 = 1,nparent
	iparent_type = my_node(icurrent_id)%parent_type(i1)
	iparent_id   = my_node(icurrent_id)%parent_id(i1)

! Prepare the inflow sets

	iadd = 0

	if ((iparent_type.eq.1).or.(iparent_type.eq.3)) then 
		call add_uncontrolled_flows(iparent_type,iparent_id,decision_var,nparam)
	end if

	! Spill is based on inflow, therfore exists for each ensemble.
	! This is currently not handled properly. 

	if(iparent_type.eq.5) then
		do j1 = 1,nwatershed
			if((parallel_track(j1)%order_type.eq.iparent_type).and. &
				(parallel_track(j1)%order_id.eq.iparent_id))call add_all_flows &
				(j1,iparent_type,iparent_id,decision_var,nparam)
		end do
		if((iparent_type.ne.iprev_type).and.(iparent_id.ne.iprev_id)) then 
			call add_controlled_flows(iparent_type,iparent_id, decision_var,nparam)
		end if
	end if

	if((iparent_type.eq.4).or.(iparent_type.eq.13)) then
		call add_controlled_flows(iparent_type,iparent_id, decision_var,nparam)
	end if
end do
15 format(a,2x, <ntime>(2x,F8.3))
! Prepare the outflow sets from uses\
! output_string = "Inflow"
! if (ifinal.eq.1) then
! 	write(35,15) output_string, (my_flow_set(iflow_set)%controlled_flows(j), j=1,ntime)
! end if
call array_ini(ntime,release, 0.0d0)
call array_ini(ntime,user_release, 0.0d0)

do i=1,nensem
	do j=1,ntime
		remaining_flow(j, i) = 0.0
	end do
end do
! call array_ini_two(nensem, ntime,remaining_flow, 0.0d0)

nchild = my_node(icurrent_id)%nchild

! Loop for adding releases to child nodes
do i1 = 1,nchild
	ichild_type = my_node(icurrent_id)%child_type(i1)
	ichild_id   = my_node(icurrent_id)%child_id(i1)
	if(ichild_type.eq.4)then
		loss_factor = my_user(ichild_id)%loss_factor
		do j  = 1,ntime
			if(ichild_type.eq.4) then 
				release (j) = release(j) + decision_var(((ichild_id-1)*ntime)+j)
				user_release(j) = decision_var(((ichild_id-1)*ntime)+j)
			end if
		end do
		! if (ifinal.eq.1) then
		! 	write(35,15) my_user(ichild_id)%name, (user_release(j), j=1,ntime)
		! end if
	end if
end do

do i = 1, nensem
	do j = 1, ntime
		remaining_flow(j, i) = my_flow_set(iflow_set)%controlled_flows(j) - release(j) &
								+ my_flow_set(iflow_set)%uncontrolled_flows(j,i)							
	end do
end do

do i1=1, nchild
	ichild_type = my_node(icurrent_id)%child_type(i1)
	ichild_id   = my_node(icurrent_id)%child_id(i1)
	if (ichild_type.eq.12) then
		if (ifinal.eq.1) then
			do i = 1, nensem
				write(35,15) my_sink(ichild_id)%name, (remaining_flow(j,i), j=1,ntime)
			end do
		end if
		! do j = 1, ntime
		! 	my_flow_set(iflow_set)%controlled_flows(j) = my_flow_set(iflow_set)%controlled_flows(j) - remaining_flow(j) - release(j)
		! end do
	end if
end do

cons_loc = nres + nuser + nres_level
do j = 1, ntime
	cons_global(cons_loc + j) = -1 !-remaining_flow(j)
end do
! Loop for flow mass balance
do j= 1,ntime
	temp = my_flow_set(iflow_set)%controlled_flows(j)
	temp = temp - release(j)
	if(temp.le.0)then 	
		my_flow_set(iflow_set)%controlled_flows(j)= 0.0
	else
		my_flow_set(iflow_set)%controlled_flows(j)= temp
	end if
end do

return
if (ifinal.eq.1) print*, 'Exiting node simul module .................'
end

subroutine add_uncontrolled_flows(iblock_type,iblock_id, decision_var,nparam)
use My_variables
implicit doubleprecision(a-h,o-z)
common /final_print/ifinal
! Loop for adding natural flow 
if (iblock_type.eq.1) then 
	do i = 1,nensem
		do j = 1,ntime
			temp = my_flow_set(iflow_set)%uncontrolled_flows(j,i) 
			temp = temp + my_watershed(iblock_id)%natural_inflows(j,i)
			my_flow_set(iflow_set)%uncontrolled_flows(j,i) = temp
		end do
	end do
else if (iblock_type.eq.3) then
	do i = 1,nensem
		do j = 1,ntime
			temp = my_flow_set(iflow_set)%uncontrolled_flows(j,i) 
			temp = temp + my_reservoir(iblock_id)%spill(j,i)
			my_flow_set(iflow_set)%uncontrolled_flows(j,i) = temp
		end do
	end do
end if

return
if (ifinal.eq.1) print *, "Exiting uncontrolled flow.........."
end


subroutine add_controlled_flows(iblock_type,iblock_id,decision_var,nparam)
use My_variables
implicit doubleprecision(a-h,o-z)
common /final_print/ifinal
double precision decision_var(nparam), loss_factor
    
! Loop for adding controlled flows
do j = 1,ntime
  	if(iblock_type.eq.4) then
		nlags = my_user(iblock_id)%nlags
		loss_factor = my_user(iblock_id)%loss_factor
    	if(j.gt.nlags)then   
			temp1 = 0.0d0
			if (nlags == 0) then
				temp1 = decision_var((iblock_id-1)*ntime+j)*(1-loss_factor)
			else 
				do k = 1,my_user(iblock_id)%nlags
					fract = (1-my_user(iblock_id)%ffraction(k))
					rflow = decision_var((iblock_id-1)*ntime+j)*fract*(1-loss_factor)
					temp1 = temp1 + rflow
				end do
			end if 
			temp = temp1
		else
	  	temp = 0.0
    	end if 
  	end if 

   if(iblock_type.eq.13) temp = my_interbasin(iblock_id)%average_flow(j)
   if(iblock_type.eq.3)  temp = decision_var((iblock_id-1)*ntime+j)
   if(iblock_type.eq.5)  temp = my_flow_set(iflow_set)%controlled_flows(j)

15	my_flow_set(iflow_set)%controlled_flows(j) = temp + my_flow_set(iflow_set)%controlled_flows(j)
end do

return
if (ifinal.eq.1) print *, 'Exiting controlled flow.......'
end

subroutine calculate_outflows(iblock_type,iblock_id,icurrent_id, icurrent_type,release,decision_var,nparam)
use My_variables
implicit doubleprecision(a-h,o-z)
common /final_print/ifinal
double precision decision_var(nparam),release(ntime)

! Loop for release calculation
	do j = 1,ntime

!       if(iblock_type.eq.4)temp = my_user(iblock_id)%demand_fract(j) &
!					*decision_var(iblock_id)
	! 	if(iblock_type.eq.3)temp = decision_var(icurrent_id+nuser*ntime)/12
	! 	if(iblock_type.eq.4)temp = decision_var((iblock_id-1)*ntime + j)
	!   if(iblock_type.eq.5)temp = decision_var(icurrent_id+nuser*ntime)/12
	! 	! if(iblock_type.eq.12)temp = decision_var(icurrent_id+nuser)/12
	! 	if(iblock_type.eq.12)temp = 0.0
		if(iblock_type.eq.4) then 
			temp = decision_var((iblock_id-1)*ntime + j)
		else
			temp = 0.0
		end if
		
   		temp1 = release(j)
		temp = temp + temp1
		release(j) = temp
	end do

return

end

! Add controlled and uncontrolled flows from another flowset
subroutine add_all_flows(iadd_set,iparent_type,iparent_id,decision_var,nparam)
use My_variables
implicit doubleprecision(a-h,o-z)

double precision decision_var(nparam)

! Add uncontrolled_flows

do j = 1,ntime
	do k = 1,nensem
		temp = my_flow_set(iadd_set)%uncontrolled_flows(j,k)
		temp1 = my_flow_set(iflow_set)%uncontrolled_flows(j,k)
		my_flow_set(iflow_set)%uncontrolled_flows(j,k) = temp + temp1
	end do
	temp = my_flow_set(iadd_set)%controlled_flows(j)
	temp1 = my_flow_set(iflow_set)%controlled_flows(j)
	my_flow_set(iflow_set)%controlled_flows(j) = temp + temp1
end do
return
end


! Check to ensure decision_var and temp_decision_var are the same
Integer function icheck_var_status(nparam, decision_var,temp_decision_var)
implicit doubleprecision (a-h,o-z)
common /final_print/ifinal
double precision decision_var(nparam),temp_decision_var(nparam)

icheck_var_status = 0
do i = 1, nparam
	if(decision_var(i).ne.temp_decision_var(i))goto 12
end do
return

12   icheck_var_status = 1

do i = 1,nparam
	temp_decision_var(i) = decision_var(i)
end do
if (ifinal.eq.1) print *,'exiting icheck_var_status..........'
return

end 
! subroutine for reservoir simulation
subroutine reser_simul(icurrent_id,q, ntime,release,storage_max, storage_ini,rate_area, &
					 simul_stor,simul_spill, simul_evapo, simul_deficit,iflag)
implicit doubleprecision(a-h,o-z)
common /final_print/ifinal
double precision q(ntime), rate_area(ntime),release(ntime) 
double precision simul_stor(ntime),simul_evapo(ntime)
double precision simul_deficit(ntime), simul_spill(ntime)

common/et_est/et_rate,storage_pre,current_flow,current_release,bal_net, iwrite

iflag = 0
storage_pre = storage_ini
! Loop for reservoir mass balance
do j = 1, ntime                  
	sum_release = release(j)                
	bal_net = q(j) - sum_release 
	current_flow  = q(j)
	current_release = sum_release
	et_rate  = rate_area(j)
	call evaporation_iter(icurrent_id, storage_current, storage_max,evapo_current, spill, deficit)
	simul_stor(j) = storage_current
	simul_evapo(j) = evapo_current
	simul_spill(j) = spill
	simul_deficit(j) = deficit
	if(deficit.gt.0.0)iflag = 1
	storage_pre = storage_current

end do
return
end


subroutine deficit_splitter(total_deficit,decision_var,nparam)
use My_variables
implicit doubleprecision(a-h,o-z)
common /final_print/ifinal
double precision total_deficit(nensem),simul_def_user(nuser), dec_var_sum
double precision deficit_split_user(nres_level,nuser),decision_var(nparam)

Integer ilevel_fail(nres_level),idef_user(nuser)

! Loop for distributing deficit among users
do k1 = 1,nensem
	call array_ini(nuser,simul_def_user,0.0d0)
	call array_ini_two(nres_level,nuser,deficit_split_user,0.0d0)
	call array_int_ini(nres_level,ilevel_fail,0)
	call array_int_ini(nuser,idef_user,0)
	if(total_deficit(k1).ne.0)then
		account_def = 0.0
		do i = 1,nres_level
			current_def = 0.0
			ilevel_fail(i) = ilevel_fail(i) + 1.0				
			do j = 1,nuser
				dec_var_sum = 0.0
				do k=1,ntime
					dec_var_sum = dec_var_sum + decision_var((j-1)*ntime + k)
				end do
				temp = my_user(j)%restr_fract(i)*dec_var_sum
				simul_def_user(j) = simul_def_user(j)+temp
				current_def = current_def+ temp			
				deficit_split_user(i,j) = temp
			end do
						
			account_def = account_def+current_def

			if(account_def.ge.total_deficit(k1))then
				ilevel = i
				adjust_def = account_def - total_deficit(k1)
				go to 17
			end if
		end do 

		go to 25

		!distribute deficits to users
		17 	do j = 1,nuser
			dec_var_sum = 0.0
			do k=1, ntime
				dec_var_sum = dec_var_sum + decision_var((j-1)*ntime + k)
			end do
			temp1 = adjust_def/current_def
			temp = my_user(j)%restr_fract(ilevel)*dec_var_sum*temp1
			simul_def_user(j) = simul_def_user(j) - temp
			deficit_split_user(ilevel,j) = my_user(j)%restr_fract(ilevel)*dec_var_sum - temp
		end do

		25 continue

		do j = 1,nuser
			if(simul_def_user(j).ge.my_user(j)%con_res_vol)idef_user(j) = idef_user(j)+ 1.0
		end do
	end if
	call get_net_ben(decision_var,nparam,idef_user,ilevel_fail,deficit_split_user,ben_net,simul_def_user)
	value_net(k1) = ben_net
end do

k2 = nres
ensem = nensem

! Loop for calculating failure probability constraints
do i = 1,nuser
    k2 = k2 +1
	if (nensem.eq.1) then
		cons_global(k2) = simul_def_user(i)
	else
		cons_global(k2) = (idef_user(i)/ensem) - my_user(i)%failure_prob
	end if
end do

! Loop for calculating target restriction constraints
do i = 1,nres_level
	k2 = k2 +1
    cons_global(k2) = (ilevel_fail(i)/ensem)-my_reservoir(1)%tar_restr_prob(i)
end do

if (ifinal.eq.1) print *, 'Exiting deficit splitter........ '
return
end

! subroutine for calculating net benefits
subroutine get_net_ben(decision_var,nparam,idef_user,ilevel_fail,deficit_split_user,ben_net,simul_def_user)
use my_variables
implicit doubleprecision(a-h,o-z)
common /final_print/ifinal
double precision decision_var(nparam),deficit_split_user(nres_level,nuser)
double precision simul_def_user (nuser)
Integer idef_user(nuser),ilevel_fail(nres_level)

ben_net = 0.0
do i = 1,nuser
	if(idef_user(i).ne.1)then
		ben_net = ben_net + my_user(i)%tariff*decision_var(i)			
	else
		temp1 = simul_def_user(i) - my_user(i)%con_res_vol
		temp = my_user(i)%penalty + temp1*my_user(i)%penalty_compen
		ben_net = ben_net - temp
  	end if
	do k = 1,nres_level
		if(ilevel_fail(k).eq.1)then
			temp = my_user(i)%res_compensation(k)*deficit_split_user(k,i)
			ben_net = ben_net - temp
		endif
	end do
end do
return
end

subroutine max_release (nparam, j, decision_var, fj)
use my_variables
implicit doubleprecision(a-h, o-z)
common /final_print/ifinal
doubleprecision decision_var(nparam), fj, x1a, x1b, x1c, x1cs
integer i, k
fj = 0.0
call constr_res(nparam,j,decision_var,gcons)
do i = 1, nuser
	do k = 1, ntime
		fj = fj + decision_var(ntime * (i - 1) + k)
	end do
end do
fj = -fj
end subroutine max_release

subroutine simple_benefits (nparam, j, decision_var, fj)
use my_variables
implicit doubleprecision(a-h, o-z)
common /final_print/ifinal
doubleprecision decision_var(nparam), fj, x1a, x1b, x1c, x1cs
integer i, k
fj = 0.0
call constr_res(nparam,j,decision_var,gcons)
do i = 1, nuser
	do k = 1, ntime
		fj = fj + my_user(i)%tariff * decision_var(ntime * (i - 1) + k)
	end do
end do
fj = -fj
end subroutine simple_benefits

subroutine spill_objective (nparam, j, decision_var, fj)
use my_variables
implicit doubleprecision(a-h, o-z)
common /final_print/ifinal
doubleprecision decision_var(nparam), fj, x1a, x1b, x1c, x1cs
integer i, k
fj = 0.0
call constr_res(nparam,j,decision_var,gcons)
do i = 1, nparam
	fj = fj + spill_values(i)
end do
end subroutine spill_objective

subroutine expected_benefits(nparam,j,decision_var,fj)
use my_variables
implicit doubleprecision(a-h,o-z)
common /final_print/ifinal
doubleprecision decision_var(nparam),fj,x1a,x1b,x1c,x1cs

fj=0.0

call constr_res(nparam,1,decision_var,gcons)

total_value = 0

do i =1, nensem
    total_value =  total_value+value_net(i)
end do

fj = -total_value

return
if (ifinal.eq.1) print *, 'Exiting expected benefits.......'
end


subroutine evaporation_iter(icurrent_id,storage_current,storage_max,evapo_current, spill,deficit) 			 
use My_variables
implicit doubleprecision(a-h,o-z)
external et_function
common/et_est/et_rate,storage_pre,current_flow,current_release,bal_net, iwrite

alpha = my_reservoir(icurrent_id)%storage_area_coeff(1)
beta = my_reservoir(icurrent_id)%storage_area_coeff(2)
 
if (storage_pre.lt.0.001) storage_pre=0.0

temp = current_flow - current_release + storage_pre
storage_min = my_reservoir(icurrent_id)%storage_min
temp_min = (storage_pre + storage_min)/2.0
rmin = storage_min + et_cal(et_rate,temp_min, alpha, beta)
temp_max = (storage_pre + storage_max)/2.0
rmax = storage_max+ et_cal(et_rate,temp_max, alpha, beta)

if (rmin.gt.temp) then
	x = storage_min
else if (rmax.lt.temp) then
	x = storage_max
else ! calculate the current storage using secant-method for root finding
	storage_current = storage_pre*0.6
	errel = 0.000001
	nvar = 1
	itmax = 5000
	icrit = 0
	x1=storage_current
	x2 = storage_pre*1.1+2 
	xacc=1.0E-06
	x=rtsec(x1,x2,xacc, alpha, beta)
end if 

! Cases for when the storage is greater than the maximum, less than the minimum
! or in between, respectively.
if (x.ge.storage_max)then
	storage_current = storage_max
	temp = 0.5*(storage_current+storage_pre)
	evapo_current = et_cal(et_rate,temp, alpha, beta)
	spill = current_flow - current_release + storage_pre - storage_max - evapo_current
	deficit = 0.0
end if 

if(x.le.storage_min)then 
	storage_current = storage_min
	temp = 0.5*(storage_current+storage_pre)
	evapo_current = et_cal(et_rate,temp, alpha, beta)
	deficit = current_release - current_flow  - storage_pre + evapo_current + storage_min
	spill = 0.0
end if 

if((x.lt.storage_max).and.(x.gt.storage_min))then
	storage_current = x
	temp = 0.5*(storage_current+storage_pre)
	evapo_current = et_cal(et_rate,temp, alpha, beta)
	deficit = 0.0
	spill = 0.0
end if 

! check mass balance
check1 = current_flow + storage_pre + deficit - spill - &
		 current_release - evapo_current - storage_current

return
end

subroutine evap_solve(storage_pre,et_rate,alpha,current_flow,current_release,x)
implicit doubleprecision(a-h,o-z)

temp1 = et_rate*alpha/2.0
temp = current_flow-current_release
temp2 = (storage_pre*(1-temp1) + temp)/(1+temp1)
x = temp2

return
end

subroutine et_function(storage_current,fnorm,nvar, alpha, beta)
implicit doubleprecision(a-h,o-z)
common/et_est/et_rate,storage_pre,current_flow,current_release,bal_net, iwrite

temp = (storage_pre+storage_current)/2.0
tmp = current_flow  - current_release + storage_pre - storage_current 
temp1 = et_cal(et_rate,temp, alpha, beta)
fnorm = tmp - temp1

return
end


double precision function et_cal(et_rate,storage, alpha, beta)
implicit doubleprecision(a-h,o-z)

! dividing by 1000000 here to go from m^3 to hm3
et_cal = et_rate*alpha*(storage**beta)/1000000

return
end

! Secant method for root finding
double precision FUNCTION rtsec(x1,x2,xacc, alpha, beta)
implicit real*8 (a-h,o-z)
INTEGER MAXIT
PARAMETER (MAXIT=30)
INTEGER j
REAL*8 dx,f,fl,swap,xl

call et_function(x1,fl,1, alpha, beta)
call et_function(x2,f,1, alpha, beta)

if (dabs(fl).lt.dabs(f)) then
	rtsec=x1
	xl=x2
	swap=fl
	fl=f
	f=swap
else
	xl=x1
	rtsec=x2
endif

do 11 j=1,MAXIT
dx=(xl-rtsec)*f/(f-fl)
xl=rtsec
fl=f
rtsec=rtsec+dx
call et_function(rtsec,f,1, alpha, beta)

if(dabs(dx).lt.xacc.or.f.eq.0.)then
	rtsec = rtsec
	return
	end if
11    continue
end


subroutine array_ini(ntime,arr, assigned_value)
implicit doubleprecision(a-h,o-z)
double precision arr(ntime),assigned_value

do i = 1, ntime
	arr(i) = assigned_value
end do

return
end				

subroutine array_ini_two(ndim1,ndim2, arr, assigned_value)
implicit doubleprecision(a-h,o-z)

doubleprecision arr(ndim1,ndim2),assigned_value

do j = 1,ndim1
	do i = 1, ndim2
	  	arr(j,i) = assigned_value
	end do
end do	

return
end				


subroutine array_int_ini(ntime,arr,my_value)
Integer arr(ntime)

do i = 1, ntime
	arr(i) = my_value
end do

return
end				
