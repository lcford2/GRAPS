! This subroutine fines upper and lower bounds of decision variables from
! relevant data  structures

subroutine  constr_decision_extract(decision_var_lb,decision_var_ub, decision_var,nparam)
USE My_variables
implicit doubleprecision(a-h,o-z)


Real*8  decision_var_lb(nparam), decision_var_ub(nparam), decision_var(nparam)

!Do i = 1,nuser

!	decision_var_lb(i) = my_user(i)%minimum_release
!	decision_var_ub(i) = my_user(i)%maximum_release
!	decision_var(i) = my_user(i)%maximum_release
!	
!End do


k =0
Do i = 1,nuser

        Do j = 1, ntime
                k = k+1
	        decision_var_lb(k) = my_user(i)%minimum_release
	        decision_var_ub(k) = my_user(i)%maximum_release
	        decision_var(k) = my_user(i)%demand_fract(j) ! read demand frac
	
        End do
End do



!Do i = (nuser*ntime)+1,(nres+nuser*ntime)!
!
!	decision_var_lb(i) = my_dir_inflows(i-nuser*ntime)%minimum_discharge   
!	decision_var_ub(i) = my_dir_inflows(i-nuser*ntime)%maximum_discharge
!	decision_var(i)    = my_dir_inflows(i-nuser*ntime)%minimum_discharge
!
!End do

if (.TRUE.) then
open(unit=43, file = 'decisionvar_details.dat',  ACTION = 'READ', STATUS = 'OLD')

DO i = 1,nparam

    read(43,*)t1
decision_var(i) = t1
END DO
close(43)
end if 

Return

end


subroutine solution_path()
Use my_variables
implicit doubleprecision (a-h,o-z)

ijump = 0

! This routine finds the solution path and order for solving the given network.
! It basically starts with a natural basin and then proceeds the order of simulation
! based on network connectivity. The simulation/mass balance is only done for 
! junction nodes and reservoirs. For rest of the system blocks, current values of the decision 
! variables are used to compute the flow through that segment.

icurrent_id = 1
icurrent_type = 1
iprev_id = 0
iprev_type = 0
iblock = 0
ijump = 0

call update_search_block(icurrent_id,icurrent_type)

DO while (ijump.eq.0)

	call child_id_finder(icurrent_id,icurrent_type)
	if(iblock.eq.ntotal_vertices)ijump =1
	icurrent_id = searched_vertices(iblock)%order_id
	icurrent_type = searched_vertices(iblock)%order_type

	
END DO 

!deallocate(searched_vertices)

return

end


subroutine child_id_finder(icurrent_id,icurrent_type)
Use my_variables
implicit doubleprecision (a-h,o-z)

	if(icurrent_type.eq.1)nchild  = my_watershed(icurrent_id)%nchild 
	if(icurrent_type.eq.3)nchild  = my_reservoir(icurrent_id)%nchild
	if(icurrent_type.eq.4)nchild  = my_user(icurrent_id)%nchild
    if(icurrent_type.eq.5)nchild  = my_node(icurrent_id)%nchild

itemp_id = icurrent_id
itemp_type = icurrent_type

DO i = 1, nchild

	icurrent_type = itemp_type
	icurrent_id = itemp_id 

	if(icurrent_type.eq.1)then
		icurrent_type = my_watershed(icurrent_id)%child_type(i)
		icurrent_id = my_watershed(icurrent_id)%child_id(i)
	else if(icurrent_type.eq.3)then
		icurrent_type = my_reservoir(icurrent_id)%child_type(i)
		icurrent_id = my_reservoir(icurrent_id)%child_id(i)
	else if(icurrent_type.eq.4)then
		icurrent_type = my_user(icurrent_id)%child_type(i)
		icurrent_id = my_user(icurrent_id)%child_id(i)
    else if(icurrent_type.eq.5)then
		icurrent_type = my_node(icurrent_id)%child_type(i)
		icurrent_id = my_node(icurrent_id)%child_id(i)
	else if(icurrent_type.eq.13)then
		icurrent_type = my_interbasin(icurrent_id)%child_type(i)
		icurrent_id = my_interbasin(icurrent_id)%child_id(i)
	end if 


	if(icurrent_type.eq.3)nparent  = my_reservoir(icurrent_id)%nparent
	if(icurrent_type.eq.4)nparent  = my_user(icurrent_id)%nparent
    if(icurrent_type.eq.5)nparent  = my_node(icurrent_id)%nparent
	if(icurrent_type.eq.12)nparent = my_sink(icurrent_id)%nparent


	if(nparent.eq.1)then

	ifound = icheck_searched(icurrent_type,icurrent_id)
	if(ifound.eq.0)then

		call update_search_block(icurrent_id,icurrent_type)
		if((icurrent_type.eq.3).or.(icurrent_type.eq.5)) &
			call update_simul_block(icurrent_id,icurrent_type)

	end if 

	else

		call subnetworks(nparent,icurrent_type,icurrent_id)

	end if

end do


return
end


recursive subroutine subnetworks(nparent,icurrent_type,icurrent_id)
Use my_variables
implicit doubleprecision(a-h,o-z)


iup_resolve = 0

Do i = 1,nparent


	if(icurrent_type.eq.3)then

		itemp_id = my_reservoir(icurrent_id)%parent_id(i)
		itemp_type = my_reservoir(icurrent_id)%parent_type(i)

     else if(icurrent_type.eq.5)then 

	 	itemp_id = my_node(icurrent_id)%parent_id(i)
		itemp_type = my_node(icurrent_id)%parent_type(i)

	else if (icurrent_type.eq.4)then

	 	itemp_id = my_user(icurrent_id)%parent_id(i)
		itemp_type = my_user(icurrent_id)%parent_type(i)

    else if (icurrent_type.eq.12) then
	 	
		itemp_id = my_sink(icurrent_id)%parent_id(i)
		itemp_type = my_sink(icurrent_id)%parent_type(i)

    end if


     if((itemp_type.eq.13).or.(itemp_type.eq.1))then
	 
	     ifound = icheck_searched(itemp_type,itemp_id)

		 if(ifound.eq.0)call update_search_block(itemp_id,itemp_type)
		 iup_resolve = iup_resolve+1


	else
	 	
		ifound = icheck_searched(itemp_type,itemp_id)

		if(ifound.eq.1)then
				iup_resolve = iup_resolve+1
		else


		if(itemp_type.eq.3)ntemp  = my_reservoir(itemp_id)%nparent
		if(itemp_type.eq.4)ntemp  = my_user(itemp_id)%nparent
		if(itemp_type.eq.5)ntemp  = my_node(itemp_id)%nparent
		if(itemp_type.eq.12)ntemp = my_sink(itemp_id)%nparent

		call subnetworks(ntemp,itemp_type,itemp_id)
		iup_resolve = iup_resolve+1

		end if

	end if

end do

if(iup_resolve.eq.nparent)then

	ifound = icheck_searched(icurrent_type,icurrent_id)
				
       if(ifound.eq.0)then

			call update_search_block(icurrent_id,icurrent_type)
			if((icurrent_type.eq.3).or.(icurrent_type.eq.5)) &
				call update_simul_block(icurrent_id,icurrent_type)

		end if

end if

return

end


		





Integer function icheck_searched(icurrent_type,icurrent_id)
Use my_variables
implicit doubleprecision (a-h,o-z)

icheck_searched = 0

DO i = 1,iblock

	if((icurrent_id.eq.searched_vertices(i)%order_id).and. &
	    (icurrent_type.eq.searched_vertices(i)%order_type))then
		
		icheck_searched =1
		GO TO 12

	end if

END DO

12	continue

return

end


subroutine update_search_block(icurrent_id,icurrent_type)
Use my_variables
implicit doubleprecision (a-h,o-z)


iblock = iblock +1


searched_vertices(iblock)%order_type = icurrent_type
searched_vertices(iblock)%order_id = icurrent_id

return
end

subroutine update_simul_block(icurrent_id,icurrent_type)
Use my_variables
implicit doubleprecision (a-h,o-z)


isimul_block = isimul_block +1


my_network_order(isimul_block)%order_type = icurrent_type
my_network_order(isimul_block)%order_id = icurrent_id

return
end
