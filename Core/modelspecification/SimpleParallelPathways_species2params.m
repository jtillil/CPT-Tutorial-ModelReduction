%%% Version: 19 Jan 2023
%%%
%%% param  = <MODELNAME>_species2params(model)
%%%
%%% This function identifies the parameters corresponding to all reactions
%%% that involve a given molecular species 
%%% 
%%% Input : model       model structure
%%%                   
%%% Output : param      strucutre parameters involved 
%%%
%%% References:
%%%
%%% 
%%% Author: Jane Knoechel and Wilhelm Huisinga
%%%

function param = SimpleParallelPathways_species2params(model)

%%% assign model indexing
I  = model.I;

param.states2Ipar = {};
param.states2nmpar = {};

% number of reactions
param.nreact = 6; % this needs to be manually supplied (not relevant here) 

%%% preprocess ODE file by replacing = by =' and ; by '; to make reactions
%%% searchable

param.states2Ipar{I.A}   = [I.k_ab I.k_ac];
param.states2Ipar{I.S}   = [I.k_ab I.k_b I.k_ac I.k_c];
param.states2Ipar{I.B}   = [I.k_ab I.k_b I.k_bd];
param.states2Ipar{I.C}   = [I.k_ac I.k_c I.k_cd];
param.states2Ipar{I.D}   = [I.k_bd I.k_cd I.k_d];

for k = 1:I.nstates
    states2param = {};
    for p = param.states2Ipar{k}
        states2param = [states2param, I.nmpar(p)];
    end
    param.states2nmpar{k} = states2param;
end

end

