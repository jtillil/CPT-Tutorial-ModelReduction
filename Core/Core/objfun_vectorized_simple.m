function vector_objvals = objfun_vectorized_simple(model, config, errtype, variability, var_obj_prctile)

[vector_objvals, ~] = objfun_vectorized(model, config, [], 120, errtype, 1, var_obj_prctile, variability.X0_pop, variability.par_pop, variability.X_ref_pop);

end
