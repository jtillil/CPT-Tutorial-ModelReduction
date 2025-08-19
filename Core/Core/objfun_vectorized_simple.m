%%% Version: October 5th, 2023
%%%
%%% call by: vector_objvals  =  redmodel_objfun_vectorized(model,mat_vector_config)
%%%
%%% This function computes errors from any reduced model. Useful for higher
%%% level algorithms. Vectorized function.
%%%
%%% Input:  model               struct specifying the model
%%%         mat_vector_config   matrix of row vectors specifying reduced model configs
%%%
%%% Output: err                 struct containing output and internal errs
%%%
%%% Citation:
%%% 
%%% ---
%%% 
%%% Authors: Johannes Tillil
%%%

function vector_objvals = objfun_vectorized_simple(model, config, errtype, variability, var_obj_prctile)

[vector_objvals, ~] = objfun_vectorized(model, config, [], 120, errtype, 1, var_obj_prctile, variability.X0_pop, variability.par_pop, variability.X_ref_pop);

end
