function state  = checkstop(params,state)
    if state.reg_obj<params.iteration
        state.stopCriterion=1;
    else
        state.stopCriterion=0;
    end
end

