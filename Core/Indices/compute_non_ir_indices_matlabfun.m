function [indx, t_ref] = compute_non_ir_indices_matlabfun(model,nmindx)

tic;
fprintf('\n Calculate %s indices supporting matlabfun\n',nmindx);

% indexing
I = model.I;

if ~ismember(nmindx,I.nondynstateclasses)
    fprintf('\n --> unsupported index type---PLEASE FIX! \n\n');
    return;
end

% check, whether classification of states is supported
if ~isempty(setxor(I.dyn,1:I.nstates))
    fprintf('\n\n --> %s indices are only computed for model with all states dynamic; PLEASE FIX! \n\n',nmindx); beep; 
    return;
end

% time vectors
t_ref  = model.t_ref; X_ref  = model.X_ref; %!% T_end = t_ref(end); 

ntstar = length(t_ref); % number (n) of tstar values

% Initialise variables necessary for analysis and keep what is already
% there
% if isfield(model,nmindx)
%     indx = model.(nmindx);
% end
% indx.index  = NaN(ntstar,I.nstates);
indx.nindex = NaN(ntstar,I.nstates);
indx.relstateerr = NaN(ntstar,I.nstates);
% tmpindex = NaN(ntstar,I.nstates);
tmpnindex = NaN(ntstar,I.nstates);
tmprelstateerr = NaN(ntstar,I.nstates);

%!% determine normalizing constant for nindex
%!% nindex_normalization = NaN(1,ntstar);
%!% for ts = 1:ntstar-1
%!%     nindex_normalization(ts) = sqrt( 1/T_end * trapz(t_ref(ts:end), X_ref(ts:end,I.output).^2) );
%!% end

% give information about progress of computation
inform = true;
if inform, fprintf(' Progress report: '); end

% determine indices of all states that are degradation products and do
% not determine indices for these states neither for input & output
% I_deg_states = find(contains(I.nmstate,'deg'));
% relevantstates = sort(setdiff(1:I.nstates,[I_deg_states I.input I.output]),'descend');
relevantstates = flip(1:I.nstates);

%%% only for testing purposes
% if model.quicktest
%     relevantstates = I.output; ntstar = 3; beep;
%     fprintf('\n --> quick test running <-- \n')
% end

% If some of the tentatively reduced models results in numerical problem,
% e.g., intergration tolerance could not be meet etc, then such a reduced
% model is not accepted. In such a case, the approximation error could
% anyway not be computed (due to the numerical issue)
o = I.output; 
for k = relevantstates
    
    any_modified_state_NaN = false;

    % set k-th species to the whichindex type and adapt the dynamic states
    model.I.(nmindx) = k; model.I.dyn = setdiff(1:I.nstates,k);
    if inform, fprintf('%d,',k); end
    
    parfor ts = 1:ntstar-1
        
        % uncomment only for testing purposes
        %if inform, fprintf('%d-',ts); end

        tstarspan = t_ref(ts:end);
        X0 = X_ref(ts,:);
        [t_per,X_per] = simModel(tstarspan, X0, model.par, model.I, model.param, model.multiple, model.odefun, model.jacfun);
        % [t_per,X_per] = model.simODE(tstarspan, X0, model.par, model);

        if any(isnan(X_per), "all")
            any_modified_state_NaN = true;
        end
        
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

            tmpnindex(ts,k)      = model.relerrnorm( tstarspan, X_per(:,o)- X_ref(ts:end,o), X_ref(ts:end,o) );
            tmprelstateerr(ts,k) = model.relerrnorm( tstarspan, X_per(:,k)- X_ref(ts:end,k), X_ref(ts:end,k) );
        else
            % some problem with integration, so just leave the initial NaN
            % values assigned
        end
    end

    if any_modified_state_NaN
        tmpnindex(:,k) = NaN;
    end
    if all(X_ref(:, k) == 0)
        tmprelstateerr(:,k) = Inf;
    end
end
indx.nindex = tmpnindex;
indx.relstateerr = tmprelstateerr;
elapsedtime = toc; fprintf('\n [elapsed time = %.1f]\n\n',elapsedtime);

%%% assign plotting routines
model.(nmindx) = indx;

% if saveresults
%     where2save = [model.savenameroot '_' nmindx '_index.mat'];
%     fprintf('  \n results saved in %s \n',where2save);
%     save(where2save,'indx','model')
% end

end
