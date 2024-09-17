%%% Version: 19 Jan 2023
%%%
%%% param = <MODELNAME>_species2params(model)
%%%
%%% This function identifies the parameters corresponding to all reactions
%%% that involve a given molecular species 
%%% 
%%% Input : model       model structure
%%%                   
%%% Output : param      strucutre parameters involved 
%%%
%%%
%%% Authors: Jane Knoechel and Wilhelm Huisinga
%%%

function param = MMEnzymeKinetics_species2params(model)

%%% assign model indexing
I  = model.I;

param.states2Ipar = {};
param.states2nmpar = {};

% number of reactions
param.nreact = 5; % this needs to be manually supplied (not relevant here) 

%%% preprocess ODE file by replacing = by =' and ; by '; to make reactions
%%% searchable

param.states2Ipar{I.A} = [I.ktrans];
param.states2Ipar{I.S}  = [I.ktrans I.kon I.koff];
param.states2Ipar{I.E}  = [I.kon I.koff I.kcat];
param.states2Ipar{I.C}  = [I.kon I.koff I.kcat];
param.states2Ipar{I.P}  = [I.kcat I.kdeg];

for k = 1:I.nstates
    states2param = {};
    for p = param.states2Ipar{k}
        states2param = [states2param, I.nmpar(p)];
    end
    param.states2nmpar{k} = states2param;
end

end

