%%% Version: 19 Jan 2023
%%%
%%% I  =  <MODELNAME>_indexing
%%%
%%% This function defines the indexing for the state variables 
%%% and parameter vector. As a result, specific state variables
%%% can be indexed in a human readable format, e.g., to
%%% access A, use X(I.A). Analogously, the parameter k_d can be 
%%% accessed via par(I.k_d) etc. 
%%%
%%% Input: none
%%%
%%% Output: I - indexing vector
%%%
%%%
%%% Authors: Jane Knoechel and Wilhelm Huisinga
%%%

function I  =  MMEnzymeKinetics_indexing

%%% define state index structure 
%%%
I.nmstate = {'A','S','E','C','P'}; % name (nm) of states
I.nstates  = length(I.nmstate);  % number (n) of states

for i = 1:I.nstates
    I.(I.nmstate{i}) = i;         % name to index mapping
end

%%% define parameter index structure
%%%
I.nmpar = {'ktrans','kon','koff','kcat','kdeg'}; % name (nm) of parameters
I.npar  = length(I.nmpar);        % number (n) of parameters

for i = 1:I.npar
        I.(I.nmpar{i}) = i;        % name to index mapping
end

%%% check for double naming and overlap between state and parameter names
%%%
if I.nstates ~= length(unique(I.nmstate))
    fprintf('--> There are double named state names---PLEASE FIX!\n\n');
    error(':-)');
end
if I.npar ~= length(unique(I.nmpar))
    fprintf('--> There are double named parameter names---PLEASE FIX!\n\n');
    error(':-)');
end
if ~isempty(intersect(I.nmstate,I.nmpar))
    fprintf('--> There is an overlap betwenn state and parameter names---PLEASE FIX!\n\n');
    error(':-)');
end

end

