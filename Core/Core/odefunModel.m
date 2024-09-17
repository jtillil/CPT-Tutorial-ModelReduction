function dX = odefunModel(X,par,I,odefun)

%%% (1) pre-processing for index analysis
%%% set negliglibe states variables to zero and account for intervention
%%%
X([I.pneg I.cneg]) = 0;

% back calculated value of states whose ODE was replaced
for p = 1:length(I.replaceODE)

    % if a state (say the k-th one) is part of I.replaceODE, then its
    % initial value X0(k) and its ODE dX(k) are replaced by the sum of
    % inital values and ODEs of the states that are specified in
    % I.replaceODEby

    % Here, this step is 'undone' by determining the value of the k-th 
    % state from the values of the states in replaceODEby
    % Enforce that difference is non-negative (ODE solver accurracy
    % might otherwise result in negative values)

    k = I.replaceODE(p);    % index of state variable

    remstates = setdiff(I.replaceODEby{p},k);
    X(k) = max(0, X(k) - sum(X(remstates)) );

end

%%% (2) call model ode
dX = odefun(X,par);

%%% (3) post-processing for index analysis
%%% set derivative of environmental & negligible states variables to zero
%%%
dX([I.env I.pneg I.cneg I.irenv_arith, I.irenv_geom I.average I.mode I.constant I.ssenv I.constregr]) = 0;

%%% determine the ODEs for the states in replaceODE
for p = 1:length(I.replaceODE)
    k  = I.replaceODE(p);       % index of state to be replaced
    states = I.replaceODEby{p}; % indices of states used for replacement
    dX(k) = sum(dX(states));
end

end