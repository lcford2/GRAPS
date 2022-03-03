
! Program used to initialize the data structures and populate
! input data for future model runs

subroutine initialize(init_nparam, init_index_cons, init_dec_vars, num_res, num_user, num_restr) !BIND(C,name='InitModule')
    ! USE ISO_C_BINDING
    Use My_variables
    USE path_vars
    implicit doubleprecision(a-h, o-z)
    common /final_print/ifinal

    double precision, allocatable :: decision_var_lb(:),decision_var_ub(:), decision_var(:)!,w(:),g(:),f(:)
    ! Integer, allocatable :: iw(:)
    integer :: init_nparam, init_index_cons, num_res, num_user, num_restr
    ! if (present(init_nparam)) then
    double precision, dimension(init_nparam), INTENT(OUT) :: init_dec_vars
    ! else
    ! double precision, allocatable :: init_dec_vars(:)
    ! end if

    ! subroutine declarations
    ! the `path.dat` file MUST be in the working directory of the caller ( either the 
    ! graps executable or the program calling the graps.so shared library )
    open(unit=9 , file = 'path.dat',                ACTION = 'READ', STATUS = 'OLD')

    0001 format(A)
    read(9, 0001)input_path
    read(9, 0001)output_path
    close(9)

    open(unit=10,  file = trim(input_path)//'input.dat',				ACTION = 'READ', STATUS = 'OLD')
    open(unit=11,  file = trim(input_path)//'watershed_details.dat',	ACTION = 'READ', STATUS = 'OLD')
    open(unit=12,  file = trim(input_path)//'nflow_details.dat',		ACTION = 'READ', STATUS = 'OLD')
    open(unit=13,  file = trim(input_path)//'reservoir_details.dat',	ACTION = 'READ', STATUS = 'OLD')
    open(unit=14,  file = trim(input_path)//'user_details.dat',		    ACTION = 'READ', STATUS = 'OLD')
    open(unit=15,  file = trim(input_path)//'node_details.dat',		    ACTION = 'READ', STATUS = 'OLD')
    open(unit=16,  file = trim(input_path)//'dir_flow_details.dat',	    ACTION = 'READ', STATUS = 'OLD')
    open(unit=17,  file = trim(input_path)//'ret_flow_details.dat',	    ACTION = 'READ', STATUS = 'OLD')
    open(unit=18,  file = trim(input_path)//'diversions_details.dat',	ACTION = 'READ', STATUS = 'OLD')
    open(unit=19,  file = trim(input_path)//'spill_flow_details.dat',	ACTION = 'READ', STATUS = 'OLD')
    open(unit=20,  file = trim(input_path)//'ibasin_flow_details.dat',	ACTION = 'READ', STATUS = 'OLD')
    open(unit=21,  file = trim(input_path)//'demand_flow_details.dat',	ACTION = 'READ', STATUS = 'OLD')
    open(unit=22,  file = trim(input_path)//'sink_details.dat',		    ACTION = 'READ', STATUS = 'OLD')
    open(unit=23,  file = trim(input_path)//'interbasin_details.dat',	ACTION = 'READ', STATUS = 'OLD')
    open(unit=25,  file = trim(input_path)//'storage_flood_rule.dat',   ACTION = 'READ', STATUS = 'OLD')
    open(unit=26,  file = trim(input_path)//'storage_supply_rule.dat',  ACTION = 'READ', STATUS = 'OLD')
    open(unit=31,  file = trim(output_path)//'storage.out')
    open(unit=54,  file = trim(output_path)//'hydro.out')
    open(unit=24,  file = trim(output_path)//'release.out')
    open(unit=28,  file = trim(output_path)//'flow_sets.out')
    open(unit=29,  file = trim(output_path)//'node_flow.out')
    open(unit=30,  file = trim(output_path)//'spill.out')
    open(unit=44,  file = trim(output_path)//'deficit.out')
    open(unit=222, file = trim(output_path)//'id_name.out')    
	open(unit=100, file = trim(output_path)//'average_spill.out')
    open(unit=101, file = trim(output_path)//'average_deficit.out')
    open(unit=102, file = trim(output_path)//'spill_prob.out')
	open(unit=103, file = trim(output_path)//'deficit_prob.out')
	open(unit=104, file = trim(output_path)//'average_hydropower.out')
    open(unit=105, file = trim(output_path)//'average_strorage.out')
    open(unit=106, file = trim(output_path)//'average_release.out')
    open(unit=107, file = trim(output_path)//'mass_balance_vars.out')
    open(unit=108, file = trim(output_path)//'res_inflow_breakdown.out')
    open(unit=110, file = trim(output_path)//'opt.out')
    open(unit=112, file = trim(output_path)//'informs.out')
    open(unit=113, file = trim(output_path)//'ignored_constraints.out')
    ! open(unit=111, file = "res_bounds.out")


    ! reading input.dat 
    read(10,*)ntime,nensem
    read(10,*)nwatershed,nnatural_flow,nres,nuser,nfnode,ndir_inflows,nret_inflows,ndiversion,nspill_flow,&
                ninterbasin_flow,ndemand_release,nsink,ninterbasin

    checksum = nwatershed+nnatural_flow+nres+nuser+nfnode+ndir_inflows+nret_inflows+ndiversion+nspill_flow+&
                ninterbasin_flow+ndemand_release+nsink+ninterbasin

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

    ! ALLOCATE space for relevant data structures.

    ! Safely allocate memory, checking because sometimes this code is run multiple times without the memory being cleared
    IF(ALLOCATED(my_user))              DEALLOCATE(my_user)
    IF(ALLOCATED(my_reservoir))         DEALLOCATE(my_reservoir)
    IF(ALLOCATED(my_node))              DEALLOCATE(my_node)
    IF(ALLOCATED(my_sink))              DEALLOCATE(my_sink)
    IF(ALLOCATED(my_watershed))         DEALLOCATE(my_watershed)
    IF(ALLOCATED(my_interbasin))        DEALLOCATE(my_interbasin)
    ALLOCATE(my_user(nuser),my_reservoir(nres),my_node(nfnode),my_sink(nsink),my_watershed(nwatershed),my_interbasin(ninterbasin))

    IF(ALLOCATED(my_dir_inflows))       DEALLOCATE(my_dir_inflows)
    IF(ALLOCATED(my_ret_inflows))       DEALLOCATE(my_ret_inflows)
    IF(ALLOCATED(my_diversions))        DEALLOCATE(my_diversions)
    ALLOCATE(my_dir_inflows(ndir_inflows),my_ret_inflows(nret_inflows),my_diversions(ndiversion))

    IF(ALLOCATED(my_spill_flow))        DEALLOCATE(my_spill_flow)
    IF(ALLOCATED(my_natural_flow))      DEALLOCATE(my_natural_flow)
    IF(ALLOCATED(my_interbasin_flow))   DEALLOCATE(my_interbasin_flow)
    IF(ALLOCATED(my_demand_release))    DEALLOCATE(my_demand_release)
    ALLOCATE(my_spill_flow(nspill_flow),my_natural_flow(nnatural_flow),my_interbasin_flow(ninterbasin_flow),&
        &my_demand_release(ndemand_release))

    ! Reads the system details from individual files until the connectivity given in input.dat comes to an end.
    icount = 0
    icount_max = checksum
    ! print *, ('-', i=1,75)
    ! print *, "Reading Input Files"
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
    ! close(24)
    ! print *, ('-', i=1,75)
    ! print *, "Done Reading Input"
    ! Convert the data structure into input parameters for the optimization routine.
    ! each reservoir should be given a minimum of one user for downstream release and it should be defined in the user_details.dat
    nparam = (nuser*ntime)
    nsimul_block = nres + nfnode
    ntotal_vertices = nwatershed + nres + nuser + nfnode + ninterbasin_flow + nsink 

    nres_level = my_reservoir(1)%nres_level
    ! total number of constraints 
    ncons = nuser + nres_level + nres

    IF(ALLOCATED(decision_var_lb))      DEALLOCATE(decision_var_lb)
    IF(ALLOCATED(decision_var_ub))      DEALLOCATE(decision_var_ub)
    IF(ALLOCATED(decision_var))         DEALLOCATE(decision_var)
    ALLOCATE(decision_var_lb(nparam), decision_var_ub(nparam), decision_var(nparam))

    IF(ALLOCATED(temp_decision_var))    DEALLOCATE(temp_decision_var)
    IF(ALLOCATED(isimul_status))        DEALLOCATE(isimul_status)
    IF(ALLOCATED(my_network_order))     DEALLOCATE(my_network_order)
    ALLOCATE(temp_decision_var(nparam),isimul_status(nsimul_block),my_network_order(nsimul_block))

    IF(ALLOCATED(searched_vertices))    DEALLOCATE(searched_vertices)
    ALLOCATE(searched_vertices(ntotal_vertices))
    IF(ALLOCATED(inproc_vertices))      DEALLOCATE(inproc_vertices)
    ALLOCATE(inproc_vertices(ntotal_vertices))

    ! lukes contraints
    IF(ALLOCATED(lukes_cons))           DEALLOCATE(lukes_cons)
    IF(ALLOCATED(cons_ignore))          DEALLOCATE(cons_ignore)
    IF(ALLOCATED(hydro_benefit))        DEALLOCATE(hydro_benefit)
    lukes_ncons = nparam * 2 + nres * 2 ! + nparam ! spill, deficit, above target storage, below target storage, hydropower max
    ALLOCATE(lukes_cons(lukes_ncons), hydro_benefit(nparam))
    ALLOCATE(cons_ignore(lukes_ncons))

    do i = 1, lukes_ncons
        lukes_cons(i) = 0.0
        cons_ignore(i) = 0
    enddo
    constraint_tolerance = 1e-2
    do i = 1, nparam
        hydro_benefit(i) = 1.0
    end do

    IF(ALLOCATED(all_release))          DEALLOCATE(all_release)
    ALLOCATE(all_release(nparam))

    opt_count = 0

    ! values for export to python
    IF(ALLOCATED(id_output))            DEALLOCATE(id_output)
    IF(ALLOCATED(cons_id))              DEALLOCATE(cons_id)
    IF(ALLOCATED(res_ids_for_spdef))    DEALLOCATE(res_ids_for_spdef)
    IF(ALLOCATED(user_id))              DEALLOCATE(user_id)
    IF(ALLOCATED(value_output))         DEALLOCATE(value_output)
    IF(ALLOCATED(spill_values))         DEALLOCATE(spill_values)
    IF(ALLOCATED(deficit_values))       DEALLOCATE(deficit_values)
    IF(ALLOCATED(cons_mag))             DEALLOCATE(cons_mag)
    IF(ALLOCATED(min_rel))              DEALLOCATE(min_rel)
    IF(ALLOCATED(max_rel))              DEALLOCATE(max_rel)
    IF(ALLOCATED(constraints))          DEALLOCATE(constraints)

    ALLOCATE(id_output(nparam), cons_id(ncons-nres_level), res_ids_for_spdef(nparam))
    ALLOCATE(user_id(nuser), value_output(nparam), spill_values(nparam), deficit_values(nparam))
    ALLOCATE(cons_mag(nres), min_rel(nuser), max_rel(nuser), constraints(ncons))
    
    do i = 1, nparam
        id_output(i) = 0
        value_output(i) = 0.0
    end do
    
    ! Initialize temporary decision variable
    call array_ini(nparam,temp_decision_var,-10.0d0)

    ! determine solution path
    call solution_path()

    call constr_decision_extract(decision_var_lb,decision_var_ub, decision_var, nparam)
    ! open(unit=60,file=trim(input_path)//'runflag.dat',ACTION = 'READ', STATUS = 'OLD')
    ! read(60,*)runflag
    ! close(60)

    
    ! allocate(initial_params(nparam+2))

    
    ! DO j = 1,74
    !     print *, initial_params(j)
    ! END DO

    ! if(runflag == 1)then 
    !     ! Read parameters for FFSQP
    !     open(unit =40, file='model_para.dat',ACTION = 'READ', STATUS = 'OLD')
    !     read(40,*)nf,mode,iprint,miter
    !     read(40,*)bigbnd,eps,epseqn,udelta 
        
    !     ! model params to modify in the above file, if needed : 
    !     ! nf : number of objective fucntions 
    !     ! mode : 110 [CBA - ref. ffsqp.f ] 
    !     ! iprint : print level information 
    !     ! miter : maximum number of iteration 
    !     ! bigbnd : plus infinity 
    !     ! eps : stopping criterion that ensures a solution, the norm of the Newton direction vector is smaller than eps 
    !     ! epseqn : tolerance of the violation of nonlinear equality constraints allowed by the user at an optimal solution 
    !     ! udelta : perturbation size 
    !     close(40)
    ! end if 

    ! output for python
    num_res = nres
    num_user = nuser
    num_restr = nres_level

    IF (ALLOCATED(cons_global))     DEALLOCATE(cons_global)
    IF (ALLOCATED(value_net))       DEALLOCATE(value_net)
    ALLOCATE(cons_global(ncons),value_net(nensem))

    IF (ALLOCATED(my_flow_set))     DEALLOCATE(my_flow_set)
    IF (ALLOCATED(parallel_track))  DEALLOCATE(parallel_track)
    ALLOCATE(my_flow_set(nwatershed),parallel_track(nwatershed))


    Do i = 1,nwatershed
        IF (ALLOCATED(my_flow_set(i)%controlled_flows))     DEALLOCATE(my_flow_set(i)%controlled_flows)
        IF (ALLOCATED(my_flow_set(i)%uncontrolled_flows))   DEALLOCATE(my_flow_set(i)%uncontrolled_flows)
        ALLOCATE(my_flow_set(i)%controlled_flows(ntime), my_flow_set(i)%uncontrolled_flows(ntime,nensem))
    end do

    ifinal = 1  
    isimul_block = 0
    init_nparam = nparam
    init_index_cons = 1
    ! init_dec_vars = decision_var
    do i=1,nparam
        init_dec_vars(i) = decision_var(i)
    end do
    
    ! DO j = 1,(nparam+2)
    !     if (j.eq.1) then
    !         initial_params(j) = nparam
    !     else if (j.eq.2) then
    !         initial_params(j) = 1 ! index_cons = 1
    !     ! else if (j.eq.3) then
    !     !     initial_params(j) = ndec_var
    !     else
    !         initial_params(j) = decision_var(j-2)
    !     end if        
    ! END DO


    21  FORMAT(I3,1x,A40)
    22 FORMAT(I3,1x,A30) 
    50 FORMAT(F10.3)
    333 FORMAT(i0,A,A)
    do j = 1, nuser
        write(222,333) j, ',', my_user(j)%name
    end do
    ! Need to flush id_name.out so it can be used by COREGS
    FLUSH(222)
    ! end of init
    ! RETURN
    ! 100	STOP
end subroutine


! Following functions are for reading input files.

! subroutine to read watershed details.dat 
Subroutine read_watershed_details(my_watershed,nwatershed,ntime, nensem)
    USE Definitions
    USE path_vars
    implicit doubleprecision(a-h,o-z)
    character*50 file_name
    character*200 path

    TYPE(watershed) my_watershed(nwatershed)


    read(11,*)inum
    read(11,20)my_watershed(inum)%name
    read(11,*)my_watershed(inum)%ID, nchild,my_watershed(inum)%drainage_area

    IF(ALLOCATED(my_watershed(inum)%child_id))          DEALLOCATE(my_watershed(inum)%child_id)
    IF(ALLOCATED(my_watershed(inum)%child_type))        DEALLOCATE(my_watershed(inum)%child_type)
    IF(ALLOCATED(my_watershed(inum)%natural_inflows))   DEALLOCATE(my_watershed(inum)%natural_inflows)

    my_watershed(inum)%nchild = nchild

    ALLOCATE(my_watershed(inum)%child_id(nchild), my_watershed(inum)%child_type(nchild))
    Do i = 1,nchild        
        read(11,*)my_watershed(inum)%child_type(i),my_watershed(inum)%child_id(i)
    END DO
    
    ALLOCATE(my_watershed(inum)%natural_inflows(ntime, nensem))

    read(11,*)file_name
    path = trim(input_path)//"InflowData/"//trim(file_name)

    open(unit=40,file=trim(path))

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

    IF (ALLOCATED(my_reservoir(inum)%child_id))         DEALLOCATE(my_reservoir(inum)%child_id)
    IF (ALLOCATED(my_reservoir(inum)%child_type))       DEALLOCATE(my_reservoir(inum)%child_type)
    IF (ALLOCATED(my_reservoir(inum)%parent_id))        DEALLOCATE(my_reservoir(inum)%parent_id)
    IF (ALLOCATED(my_reservoir(inum)%parent_type))      DEALLOCATE(my_reservoir(inum)%parent_type)
    IF (ALLOCATED(my_reservoir(inum)%spill_type))       DEALLOCATE(my_reservoir(inum)%spill_type)
    IF (ALLOCATED(my_reservoir(inum)%crest_level))      DEALLOCATE(my_reservoir(inum)%crest_level)
    IF (ALLOCATED(my_reservoir(inum)%discharge_max))    DEALLOCATE(my_reservoir(inum)%discharge_max)
    IF (ALLOCATED(my_reservoir(inum)%elevation_rvalve)) DEALLOCATE(my_reservoir(inum)%elevation_rvalve)
    IF (ALLOCATED(my_reservoir(inum)%area_rvalve))      DEALLOCATE(my_reservoir(inum)%area_rvalve)
    IF (ALLOCATED(my_reservoir(inum)%loss_coeff_max))   DEALLOCATE(my_reservoir(inum)%loss_coeff_max)
    IF (ALLOCATED(my_reservoir(inum)%loss_coeff_min))   DEALLOCATE(my_reservoir(inum)%loss_coeff_min)
    IF (ALLOCATED(my_reservoir(inum)%rule_curve_upper)) DEALLOCATE(my_reservoir(inum)%rule_curve_upper)
    IF (ALLOCATED(my_reservoir(inum)%rule_curve_lower)) DEALLOCATE(my_reservoir(inum)%rule_curve_lower)
    IF (ALLOCATED(my_reservoir(inum)%evaporation_rate)) DEALLOCATE(my_reservoir(inum)%evaporation_rate)
    IF (ALLOCATED(my_reservoir(inum)%tar_restr_prob))   DEALLOCATE(my_reservoir(inum)%tar_restr_prob)

    read(13,20)my_reservoir(inum)%name 
    read(13,*)my_reservoir(inum)%latitude, my_reservoir(inum)%longitude
    ! print *, my_reservoir(inum)%name
    read(13,*)my_reservoir(inum)%elevation_max, my_reservoir(inum)%elevation_min 
    read(13,*)my_reservoir(inum)%storage_max, my_reservoir(inum)%storage_min, my_reservoir(inum)%current_storage
    ! write(111, "(A, F, F)") my_reservoir(inum)%name, my_reservoir(inum)%storage_max, my_reservoir(inum)%storage_min
    read(13,*)my_reservoir(inum)%elevation_storage_coeff(1),my_reservoir(inum)%elevation_storage_coeff(2) &
            ,my_reservoir(inum)%elevation_storage_coeff(3)
    read(13,*)my_reservoir(inum)%storage_area_coeff(1),my_reservoir(inum)%storage_area_coeff(2)

    read(13,*)my_reservoir(inum)%nspillway,my_reservoir(inum)%number_outlets
    read(13,*)my_reservoir(inum)%nres_level
    read(13,*)my_reservoir(inum)%nchild,my_reservoir(inum)%nparent

    n1 = my_reservoir(inum)%nchild
    n2 = my_reservoir(inum)%nparent
    n4 = my_reservoir(inum)%nspillway

    ALLOCATE(my_reservoir(inum)%child_id(n1),my_reservoir(inum)%child_type(n1))
    ALLOCATE(my_reservoir(inum)%parent_id(n2),my_reservoir(inum)%parent_type(n2))
    ALLOCATE(my_reservoir(inum)%spill_type(n4),my_reservoir(inum)%crest_level(n4),my_reservoir(inum)%discharge_max(n4))

    Do i = 1,n4
        read(13,*)my_reservoir(inum)%spill_type(i),my_reservoir(inum)%crest_level(i), my_reservoir(inum)%discharge_max(i)
    End do

    Do i = 1,n1
        read(13,*)my_reservoir(inum)%child_type(i), my_reservoir(inum)%child_id(i)
    end do

    Do i = 1,n2
        read(13,*)my_reservoir(inum)%parent_type(i), my_reservoir(inum)%parent_id(i)
    end do

    n3 = my_reservoir(inum)%number_outlets

    ALLOCATE(my_reservoir(inum)%elevation_rvalve(n3),my_reservoir(inum)%area_rvalve(n3))
    ALLOCATE(my_reservoir(inum)%loss_coeff_max(n3),my_reservoir(inum)%loss_coeff_min(n3))

    Do i = 1, n3
        read(13,*)t1,t2,t3,t4,t5,t6
        my_reservoir(inum)%elevation_rvalve(i) = t1
        my_reservoir(inum)%area_rvalve(i)	   = t2
        my_reservoir(inum)%loss_coeff_max(i)   = t3
        my_reservoir(inum)%loss_coeff_max(i)   = t4
        my_reservoir(inum)%target_storage      = t5
        my_reservoir(inum)%storage_prob        = t6
    End do

    
    ALLOCATE(my_reservoir(inum)%rule_curve_upper(ntime), &
             my_reservoir(inum)%rule_curve_lower(ntime), & 
             my_reservoir(inum)%evaporation_rate(ntime))

    read(13,*)(my_reservoir(inum)%rule_curve_upper(i),i=1,ntime)
    read(13,*)(my_reservoir(inum)%rule_curve_lower(i),i=1,ntime)
    read(13,*)(my_reservoir(inum)%evaporation_rate(i),i=1,ntime)

    n4 = my_reservoir(inum)%nres_level

    ALLOCATE(my_reservoir(inum)%tar_restr_prob(n4))
        
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

    IF (ALLOCATED(my_user(inum)%child_id))          DEALLOCATE(my_user(inum)%child_id)
    IF (ALLOCATED(my_user(inum)%child_type))        DEALLOCATE(my_user(inum)%child_type)
    IF (ALLOCATED(my_user(inum)%parent_id))         DEALLOCATE(my_user(inum)%parent_id)
    IF (ALLOCATED(my_user(inum)%parent_type))       DEALLOCATE(my_user(inum)%parent_type)
    IF (ALLOCATED(my_user(inum)%demand_fract))      DEALLOCATE(my_user(inum)%demand_fract)
    IF (ALLOCATED(my_user(inum)%restr_fract))       DEALLOCATE(my_user(inum)%restr_fract)
    IF (ALLOCATED(my_user(inum)%res_compensation))  DEALLOCATE(my_user(inum)%res_compensation)
    IF (ALLOCATED(my_user(inum)%tail_elevation))    DEALLOCATE(my_user(inum)%tail_elevation)
    IF (ALLOCATED(my_user(inum)%ffraction))         DEALLOCATE(my_user(inum)%ffraction)

    read(14,20)my_user(inum)%name
    read(14,*)my_user(inum)%id, my_user(inum)%user_type,my_user(inum)%nchild, my_user(inum)%nparent,my_user(inum)%nres_level

    n1 = my_user(inum)%nchild
    n2 = my_user(inum)%nparent
    nres_level = my_user(inum)%nres_level

    ALLOCATE(my_user(inum)%child_id(n1),my_user(inum)%child_type(n1))
    ALLOCATE(my_user(inum)%parent_id(n2),my_user(inum)%parent_type(n2))
    ALLOCATE(my_user(inum)%demand_fract(ntime), my_user(inum)%restr_fract(nres_level), my_user(inum)%res_compensation(nres_level))

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


        ALLOCATE(my_user(inum)%tail_elevation(ntime))
        ! Don't really like this changing every ensemble
        ! I don't know if it has to change but I think it would make the user_details file
        ! extremely long when using ensembles
        DO k = 1,nensem
            read(14,*)(my_user(inum)%tail_elevation(j1), j1=1,ntime)
        END DO
    end if

    read(14,*)nlags

    my_user(inum)%nlags = nlags

    ALLOCATE(my_user(inum)%ffraction(nlags))

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

    IF (ALLOCATED(my_node(inum)%child_id))    DEALLOCATE(my_node(inum)%child_id)
    IF (ALLOCATED(my_node(inum)%child_type))  DEALLOCATE(my_node(inum)%child_type)
    IF (ALLOCATED(my_node(inum)%parent_id))   DEALLOCATE(my_node(inum)%parent_id)
    IF (ALLOCATED(my_node(inum)%parent_type)) DEALLOCATE(my_node(inum)%parent_type)


    read(15,20)my_node(inum)%name

    read(15,*)my_node(inum)%id, my_node(inum)%nchild, my_node(inum)%nparent

    n1 = my_node(inum)%nchild
    n2 = my_node(inum)%nparent


    ALLOCATE(my_node(inum)%child_id(n1),my_node(inum)%child_type(n1))
    ALLOCATE(my_node(inum)%parent_id(n2),my_node(inum)%parent_type(n2))

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

    IF (ALLOCATED(my_sink(inum)%parent_id))         DEALLOCATE(my_sink(inum)%parent_id)
    IF (ALLOCATED(my_sink(inum)%parent_type))       DEALLOCATE(my_sink(inum)%parent_type)


    read(22,*)my_sink(inum)%name

    read(22,*)my_sink(inum)%id, my_sink(inum)%nparent

    n2 = my_sink(inum)%nparent


    ALLOCATE(my_sink(inum)%parent_id(n2),my_sink(inum)%parent_type(n2))

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
    
    IF (ALLOCATED(my_interbasin(inum)%child_id))          DEALLOCATE(my_interbasin(inum)%child_id)
    IF (ALLOCATED(my_interbasin(inum)%child_type))        DEALLOCATE(my_interbasin(inum)%child_type)
    IF (ALLOCATED(my_interbasin(inum)%average_flow))        DEALLOCATE(my_interbasin(inum)%average_flow)

    read(23,20)my_interbasin(inum)%name
    read(23,*) my_interbasin(inum)%drainage_area

    read(23,*)my_interbasin(inum)%id, my_interbasin(inum)%nchild

    n1 = my_interbasin(inum)%nchild


    ALLOCATE(my_interbasin(inum)%child_id(n1),my_interbasin(inum)%child_type(n1))

    Do i = 1,my_interbasin(inum)%nchild

        read(23,*)my_interbasin(inum)%child_type(i), my_interbasin(inum)%child_id(i)

    end do

    ALLOCATE(my_interbasin(inum)%average_flow(ntime))

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

! subroutine array_ini(ntime,arr, assigned_value)
! implicit doubleprecision(a-h,o-z)
! double precision arr(ntime),assigned_value

!     do i = 1, ntime
!         arr(i) = assigned_value
!     end do

! return
! end	

subroutine  constr_decision_extract(decision_var_lb,decision_var_ub, decision_var, nparam)
    USE My_variables
    USE path_vars
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
        open(unit=43, file = trim(input_path)//'decisionvar_details.dat',  ACTION = 'READ', STATUS = 'OLD')
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

! This routine finds the solution path and order for solving the given network.
! It basically starts with a natural basin and then proceeds the order of simulation
! based on network connectivity. The simulation/mass balance is only done for 
! junction nodes and reservoirs. For rest of the system blocks, current values of the decision 
! variables are used to compute the flow through that segment.

    ijump = 0
    icurrent_id = 1
    icurrent_type = 1
    iprev_id = 0
    iprev_type = 0
    iblock = 0
    ijump = 0

    inproc_top = 0
    do i = 1, ntotal_vertices
        inproc_vertices(i)%order_id = 0
        inproc_vertices(i)%order_type = 0
    end do


    call update_search_block(icurrent_id,icurrent_type)
    DO while (ijump.eq.0)
        call push_inproc_vert(icurrent_id, icurrent_type)
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
            call pop_inproc_vert()
            ifound = icheck_searched(icurrent_type,icurrent_id)
            if(ifound.eq.0)then
                call update_search_block(icurrent_id,icurrent_type)
                if((icurrent_type.eq.3).or.(icurrent_type.eq.5)) &
                    call update_simul_block(icurrent_id,icurrent_type)
            end if 
        else
            call push_inproc_vert(icurrent_id, icurrent_type)
            call subnetworks(nparent,icurrent_type,icurrent_id)
            call pop_inproc_vert()
        end if
    end do
    return  
end


recursive subroutine subnetworks(nparent,icurrent_type,icurrent_id)
Use my_variables
implicit doubleprecision(a-h,o-z)
integer i_inproc

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
            if(ifound.eq.0) then 
                call update_search_block(itemp_id,itemp_type)
            end if
            iup_resolve = iup_resolve+1
        else            
            ifound = icheck_searched(itemp_type,itemp_id)
            i_inproc = icheck_inproc(itemp_type, itemp_id)
            if(ifound.eq.1)then
                iup_resolve = iup_resolve+1
            else if (i_inproc.eq.1) then
                iup_resolve = iup_resolve + 1
            else
                if(itemp_type.eq.3)ntemp  = my_reservoir(itemp_id)%nparent
                if(itemp_type.eq.4)ntemp  = my_user(itemp_id)%nparent
                if(itemp_type.eq.5)ntemp  = my_node(itemp_id)%nparent
                if(itemp_type.eq.12)ntemp = my_sink(itemp_id)%nparent
                call push_inproc_vert(itemp_id, itemp_type)
                call subnetworks(ntemp,itemp_type,itemp_id)
                call pop_inproc_vert()
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

Integer function icheck_inproc(icurrent_type,icurrent_id)
    Use my_variables
    implicit doubleprecision (a-h,o-z)

    icheck_inproc = 0
    DO i = 1,iblock
        if((icurrent_id.eq.inproc_vertices(i)%order_id).and. &
            (icurrent_type.eq.inproc_vertices(i)%order_type))then            
            icheck_inproc =1
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

subroutine push_inproc_vert(icurrent_id,icurrent_type)
    Use my_variables
    implicit doubleprecision (a-h,o-z)

    inproc_top = inproc_top + 1
    inproc_vertices(inproc_top)%order_type = icurrent_type
    inproc_vertices(inproc_top)%order_id = icurrent_id

    return
end

subroutine pop_inproc_vert()
    Use my_variables
    implicit doubleprecision (a-h,o-z)

    inproc_vertices(inproc_top)%order_type = 0
    inproc_vertices(inproc_top)%order_id = 0
    inproc_top = inproc_top - 1

    return
end

subroutine print_inproc_vert()
    Use my_variables
    implicit doubleprecision (a-h,o-z)

    do i = 1, inproc_top
        print *, i, inproc_vertices(i)%order_type, &
                    inproc_vertices(i)%order_id
    end do

    return
end

subroutine remove_inproc_vert(icurrent_id,icurrent_type)
    Use my_variables
    implicit doubleprecision (a-h,o-z)
    integer rm_index

    rm_index = inproc_top
    do i = 1, inproc_top
        if((icurrent_id.eq.inproc_vertices(i)%order_id).and. &
            (icurrent_type.eq.inproc_vertices(i)%order_type)) then

                inproc_vertices(i)%order_type = 0
                inproc_vertices(i)%order_id = 0
                rm_index = i
            GO TO 13
        end if
    end do
    13 continue
    do i = rm_index, inproc_top - 1
        inproc_vertices(i)%order_type = inproc_vertices(i+1)%order_type
        inproc_vertices(i)%order_id = inproc_vertices(i+1)%order_id
    end do

    inproc_top = inproc_top - 1

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
