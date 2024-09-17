%%% Version: 19 Jan 2023
%%%
%%% call by: default_state2linestyle(I)
%%%
%%% This function assigns the line style to the given state variable to 
%%% facilitate plotting and to have always the same style for the specific 
%%% species.
%%%
%%% 
%%% Authors: Jane Knoechel and Wilhelm Huisinga
%%%

function state2style = default_state2linestyle(I)

%%%
usedefaultlinestyle = true;

if usedefaultlinestyle

    %%% define when colormap repeats the colors so that a different line styles
    %%% should be used then
    repeat = 7;

    styles = {'-','--','-.',':'};
    nstyles = length(styles);

    % initialize output and assign colors
    state2style = cell(I.nstates,1);

    for k = 1:I.nstates
        state2style(k) = styles(mod(floor(k/repeat),nstyles)+1);
    end

else
    %%% allows to specify line styles that facilitate reading when used 
    %%% in plots
    %%% if desired, copy file into local modelspecification folder, 
    %%% rename to <PROJECT NAME>_state2linestyle.m; change
    %%% in line 21 above usedefaultnames to FALSE, and change in the local  
    %%% <PROJECT NAME>_model_set_up_details.m file the variable to 
    %%% model.setup.linestyleprovided to TRUE

    %%% -----------------------------------------------------------------------
    % state2style = ...;
    %%% -----------------------------------------------------------------------
end
