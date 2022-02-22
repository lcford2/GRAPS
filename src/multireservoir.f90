Program multireservoir
    ! USE init 
    ! USE all_simul  
    USE path_vars  
    implicit doubleprecision(a-h, o-z)
    common /final_print/ifinal
    integer :: nparam, index_cons, num_res, num_user, num_restr, runflag, ncons, lukes_ncons
    double precision, allocatable :: decision_vars(:), value_output(:), constraints(:)
    double precision, allocatable :: cons_mag(:), min_rel(:), max_rel(:), spill_values(:), deficit_values(:)
    double precision :: gcons
    character(LEN=3) :: func_flag
    integer , allocatable :: id_output(:), cons_id(:), user_id(:), res_ids_for_spdef(:)
    
    open(unit=9 , file = 'path.dat',                ACTION = 'READ', STATUS = 'OLD')

    0001 format(A)
    read(9, 0001)input_path
    read(9, 0001)output_path
    close(9)

    open(unit=60,file=trim(input_path)//'runflag.dat',ACTION = 'READ', STATUS = 'OLD')
    read(60,*)runflag
    close(60)

    open(unit=10, file = trim(input_path)//'input.dat',				ACTION = 'READ', STATUS = 'OLD')
    read(10,*) ntime, nensem
    read(10,*)nwatershed,nnatural_flow,num_res,num_user,nfnode,ndir_inflows,nret_inflows,ndiversion,nspill_flow,ninterbasin_flow,ndemand_release,nsink,ninterbasin
    nparam = (num_user*ntime)
    close(10)

    Allocate(decision_vars(nparam), id_output(nparam), value_output(nparam))
    

    call initialize(nparam,index_cons, decision_vars, num_res, num_user, num_restr)
    ncons = num_user+num_restr+num_res
    lukes_ncons = nparam * 2 + num_res
    Allocate(constraints(ncons), cons_id(num_res+num_user), cons_mag(num_res))
    Allocate(min_rel(num_user), max_rel(num_user), user_id(num_user))
    Allocate(spill_values(nparam), deficit_values(nparam), res_ids_for_spdef(nparam))

    if (runflag.eq.0) then
        call constr_res(nparam,index_cons,decision_vars,gcons)!,id_output,value_output, &
                    ! constraints,cons_id,cons_mag,min_rel,max_rel,user_id, &
                    ! spill_values, deficit_values, res_ids_for_spdef)
    else
        func_flag = "mhb"
        call optimize(nparam, lukes_ncons, decision_vars, func_flag)
    end if
    
end program multireservoir
