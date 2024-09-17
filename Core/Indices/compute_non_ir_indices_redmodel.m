%%% Version: 19 Jan 2023
%%%
%%% call by: indx = compute_non_ir_indices_redmodel(model,whichindex)
%%%
%%% This function computes an "perturbation" index for a reduced model;
%%% it particularly allows to determine an quasi-steady state index
%%% based on a partial steady state in addition to conservation laws and
%%% transformed ODEs
%%%
%%% Input:  states              to compute the index specified by
%%%         model               structure specifying the model
%%%         saveresults         1 = yes, 0 = no
%%%
%%% Output: indx                index structure with subfield
%%%                             index, nindex,relstateerr
%%%                             and plotting functions
%%%
%%%
%%% Authors: Jane Knoechel and Wilhelm Huisinga
%%%

function indx = compute_non_ir_indices_redmodel(states,Ired,model,saveresults)

tic;
fprintf('\n Calculate index for reduced model \n');

%%% set up perturbed (reduced) model
I = model.I; sumofsizes = 0; I_nondynstates = []; 

for s = I.nondynstateclasses
    class = s{:};
    if isfield(Ired,class)
        % sum of sizes needed to check below for pairwise disjointness
        sumofsizes = sumofsizes + length(Ired.(class));
        I_nondynstates = [I_nondynstates Ired.(class)];
    else
        Ired.(class) = [];
    end   
end
Ired.dyn = setdiff(1:I.nstates,I_nondynstates);

%%% consistency check for Ired: should be pairwise disjoint; use the fact: 
%%% finite sets are pairwise disjoint if size of their union equals 
%%% the sum of their sizes
if ~( length(unique(I_nondynstates)) == sumofsizes )
    fprintf('\n\n --> Set specified in Ired are not pairwise disjoint; PLEASE FIX! \n\n'); beep; 
    return;
end

%%% enrich a Ired by all fields of I that are not yet field of Ired
for fnm = fieldnames(I)'
    if ~isfield(Ired,fnm{:})
        Ired.(fnm{:}) = I.(fnm{:});
    end
end
redmodel = model; redmodel.I = Ired;

% time vectors
t_ref  = model.t_ref; X_ref  = model.X_ref; %!% T_end = t_ref(end); 

ntstar = length(t_ref); % number (n) of tstar values

% Initialise variables necessary for analysis 
indx.nindex = NaN(ntstar,I.nstates);
indx.relstateerr = NaN(ntstar,I.nstates);

%!% determine normalizing constant for nindex
%!% nindex_normalization = NaN(1,ntstar);
%!% for ts = 1:ntstar-1
%!%     nindex_normalization(ts) = sqrt( 1/T_end * trapz(t_ref(ts:end), X_ref(ts:end,I.output).^2) );
%!% end

% give information about progress of computation
inform = true;
if inform, fprintf(' Progress report: '); end

%%% only for testing purposes
if model.quicktest
    states = I.output; ntstar = 3; beep;
    fprintf('\n --> quick test running <-- \n')
end

% If some of the tentatively reduced models results in numerical problem,
% e.g., intergration tolerance could not be meet etc, then such a reduced
% model is not accepted. In such a case, the approximation error could
% anyway not be computed (due to the numerical issue)
o = I.output; 
for k = states
    
    if inform, fprintf('%d,',k); end
    
    for ts = 1:ntstar-1
        
        % uncomment only for testing purposes
        %if inform, fprintf('%d-',ts); end

        tstarspan = t_ref(ts:end);
        X0 = X_ref(ts,:);
        [t_per,X_per] = model.simODE(tstarspan, X0, redmodel.par, redmodel);
                
        if ~isnan(t_per)           
            % numerics seems to be fine, so compute the approximation error
            % if tstarspan only contains two values, ODEs interpretes
            % it differently, so reduce output to these to values
            if length(tstarspan)==2
                X_per = X_per([1 end],:);
            end
%!%             indx.index(ts,k)  = sqrt( 1/T_end * trapz(tstarspan, ( X_per(:,I.output)- X_ref(ts:end,I.output) ).^2) );
%!%             indx.nindex(ts,k) = indx.index(ts,k)/nindex_normalization(ts);
%!%             indx.relstateerr(ts,k) = model.errfun( tstarspan, X_ref(ts:end,k), X_per(:,k) );

            indx.nindex(ts,k)      = model.relerrnorm( tstarspan, X_per(:,o)- X_ref(ts:end,o), X_ref(ts:end,o) );
            indx.relstateerr(ts,k) = model.relerrnorm( tstarspan, X_per(:,k)- X_ref(ts:end,k), X_ref(ts:end,k) );
        else
            % some problem with integration, so just leave the initial NaN
            % values assigned
        end
    end
end
elapsedtime = toc; fprintf('\n [elapsed time = %.1f]\n\n',elapsedtime);

%%% assign plotting routines
model.red = indx;

if saveresults
    where2save = [model.savenameroot '_red_index.mat'];
    fprintf('  \n results saved in %s \n',where2save);
    save(where2save,'indx','model')
end


%%% during testing only
figure(5); hold on; plot(model.t_ref,model.red.nindex(:,I.C),'k--','LineWidth',2,'DisplayName','tpss'); hold off;
figure(6); hold on; plot(model.t_ref,model.red.relstateerr(:,I.C),'k--','LineWidth',2,'DisplayName','tpss'); hold off;




end

