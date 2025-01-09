%%% Version: July 11th, 2022
%%%
%%% I  =  Wajima2009BloodCoagulation_indexing
%%%
%%% This function defines the indexing for the state variables and parameter vector of 
%%% the blood coagulation model published by Wajima et al 2009
%%%
%%% Parameters are considered as states for observability calculation
%%%
%%% Input: none
%%%
%%% Output: I - indexing vector
%%%
%%% 
%%% Authors: Undine Falkenhagen, Jane Knoechel and Wilhelm Huisinga
%%%


function I = Gulati2014BClumped_indexing

%%% define state index structure 
%%%
I.nmstate = {'AVenom','CVenom', 'lumped', 'IIa','Fg'}; %,'lump'}; %,'PT'

I.nstates  = length(I.nmstate);  % number (n) of states

for i = 1:I.nstates
    I.(I.nmstate{i}) = i;         % name to index mapping
end

%%% define parameter index structure
%%%
I.nmpar = {'dIIa', 'vlumpedVenom2IIa', 'klumpedVenom2IIa', ...
    'dFg', 'pFg', 'vFgIIa2F', 'kFgIIa2F', ...
    'dVenom', 'kaVenom', 'plumped', 'dlumped'};

I.npar  = length(I.nmpar);        % number (n) of parameters

for i = 1:I.npar
        I.(I.nmpar{i}) = i;        % name to index mapping
end


%%% check for overlap between state and parameter names
%%%
if ~isempty(intersect(I.nmstate,I.nmpar))
    fprintf('--> There is an overlap betwenn state and parameter names :-( Please fix!\n\n');
    error('--> Ending here');
end


end
