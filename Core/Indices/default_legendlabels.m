%%% Version: 19 Jan 2023
%%%
%%% I  = default_legendlabels(model)
%%%
%%% This function assigns the better readable labels to states
%%% 
%%% Input:  model       model structure
%%%                   
%%% Output: I           updated strucutre 
%%% 
%%%
%%% Author: Jane Knoechel and Wilhelm Huisinga
%%%


function I = default_legendlabels(I)

%%%
usedefaultnames = true;

S = struct;
if usedefaultnames
    %%% use the same name as used in the index structure
    for k = 1:I.nstates
        nm = I.nmstate{k};
        S.(nm) = nm;
    end
else
    %%% allows to specify different name that facilitate reading when used 
    %%% in legends
    %%% if desired, copy file into local modelspecification folder, 
    %%% rename to <PROJECT NAME>_legendlabels.m; change
    %%% in line 19 above usedefaultnames to FALSE, and change in the local  
    %%% <PROJECT NAME>_model_set_up_details.m file the variable to 
    %%% model.setup.legendlabelsprovided to TRUE

    %%% -----------------------------------------------------------------------
    S = struct(...
        'A','my name for A',...
        'B','my name for B',...
        'C','my name for C');
    %%% -----------------------------------------------------------------------
end


for k = 1:I.nstates
    I.nmstatelegend{k} = S.(I.nmstate{k});
end
