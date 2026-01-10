function model = compute_and_analyse_indices_matlabfun(model,action)

fprintf('\n\n  This is compute_and_analyse_indices ...\n\n')

switch action
    case 'compute'
        saveresults = 1;      
        [model.ir, model.contr, model.obs] = compute_ir_indices_matlabfun(model,saveresults);

        % compute env, neg, pss and int indices
        model.env  = compute_non_ir_indices(model,'env',saveresults);
        model.pneg = compute_non_ir_indices(model,'pneg',saveresults);
        model.cneg = compute_non_ir_indices(model,'cneg',saveresults);
        model.pss  = compute_non_ir_indices(model,'pss',saveresults);
        
    case 'load'
        model = load_indices(model); 
    otherwise
        fprintf('\n --> unsupported action type; only ''compute'' or ''load'' allowed---PLEASE FIX! \n\n');
        return;
end

%%% assign plotting routines
model.ir.plotnindex   = @(figNr,varargin) plot_ir_nindex(figNr,model,'ir',varargin{:});
model.contr.plotindex = @(figNr,varargin) plot_ir_index(figNr,model,'contr',varargin{:});
model.obs.plotindex   = @(figNr,varargin) plot_ir_index(figNr,model,'obs',varargin{:});

for whichindex = {'env','pneg','cneg','pss'}
    nmindx = whichindex{:};
    model.(nmindx).plotnindex      = @(figNr,varargin) plot_non_ir_nindex(figNr,model,nmindx,varargin{:});
    model.(nmindx).plotrelstateerr = @(figNr,varargin) plot_non_ir_relstateerr(figNr,model,nmindx,varargin{:});
end

model.threshold = 0.1; % threshold for nindices and relstateerr

%%% plot number of nir-indices above the threshold
model.plotnnirabovethreshold = @(figNr) plot_number_nir_above_threshold(figNr,model);

%%% plotting all indices and rel state error for a given state
model.plotnindices4state = @(figNr,I_state) plotnindices4state(figNr,I_state,model);
model.plotrelstateerr4state = @(figNr,I_state) plotrelstateerr4state(figNr,I_state,model);


% analyse the indices (incl. normalized ones)
model.analysis = analyse_all_indices(model);

end

%%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%%% LOCAL SUB-ROUTINES
%%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% -------------------------------------------------------------------------
function [] = plot_ir_index(figNr,model,whichindex,varargin)

I = model.I;
states2plot = 1:I.nstates;

if nargin == 4
    states2plot = varargin{1};
end
state2col = model.state2col;
state2linestyle = model.state2linestyle;
index = model.(whichindex).index;

%%% potential time unit conversion
unit = model.unit;
t_ref = unit.graphic.transftime(model.t_ref);

figure(figNr); clf
for k = 1:length(states2plot)
    s = states2plot(k);
    plot(t_ref,index(:,states2plot(k)),'LineWidth',2,'LineStyle',state2linestyle{s},'Color',state2col(s,:));
    hold on;
end
xlabel(sprintf('t [%s]',unit.graphic.time)); ylabel(sprintf('%s-index',whichindex)); title(sprintf('%s-index',whichindex))
set(gca,'yscale','lin'); legend(I.nmstatelegend(states2plot),'Location','eastoutside');
% zoom in/out on the Y axis by changing the command below
%ylim([1e2 5e7])
model.makeplotbold(figNr); drawnow;

end

% -------------------------------------------------------------------------
function [] = plot_ir_nindex(figNr,model,whichindex,varargin)

I = model.I;
states2plot = 1:I.nstates;

if nargin == 4
    states2plot = varargin{1};
end
state2col = model.state2col;
state2linestyle = model.state2linestyle;
nindex = model.(whichindex).nindex;

%%% potential time unit conversion
unit = model.unit;
t_ref = unit.graphic.transftime(model.t_ref);

figure(figNr); clf
for k = 1:length(states2plot)
    s = states2plot(k);
    plot(t_ref,nindex(:,states2plot(k)),'LineWidth',2,'LineStyle',state2linestyle{s},'Color',state2col(s,:));
    hold on;
end
xlabel(sprintf('t [%s]',unit.graphic.time)); ylabel(sprintf('n%s-index',whichindex)); title(sprintf('n%s-index',whichindex))
set(gca,'yscale','lin'); legend(I.nmstatelegend(states2plot),'Location','eastoutside');
% zoom in/out on the Y axis by changing the command below
xlim(model.unit.graphic.xlim); ylim([-0.01 1])
model.makeplotbold(figNr); drawnow;

end

% -------------------------------------------------------------------------
function [] = plot_non_ir_nindex(figNr,model,whichindex,varargin)

I = model.I;

%%% states for graphical output
state2col = model.state2col;
state2linestyle = model.state2linestyle;
states2plot = 1:I.nstates;
if nargin == 4
    states2plot = varargin{1};
end

%%% quantity of interest
nindex = model.(whichindex).nindex;

%%% potential time unit conversion
unit = model.unit;
t_ref = unit.graphic.transftime(model.t_ref);

%%% graphical output
figure(figNr); clf
for k = 1:length(states2plot)
    s = states2plot(k);
    plot(t_ref,nindex(:,states2plot(k)),'LineWidth',2,'LineStyle',state2linestyle{s},'Color',state2col(s,:));
    hold on;
end
xlabel(sprintf('t [%s]',unit.graphic.time)); ylabel(sprintf('%s-nindex',whichindex)); title(sprintf('%s-nindex',whichindex))
set(gca,'yscale','lin'); legend(I.nmstatelegend(states2plot),'Location','eastoutside');
xlim(model.unit.graphic.xlim); ylim([-0.01 1])
model.makeplotbold(figNr); drawnow;

end

% -------------------------------------------------------------------------
function [] = plot_non_ir_relstateerr(figNr,model,whichindex,varargin)

I = model.I;

%%% states for graphical output
state2col = model.state2col;
states2plot = 1:I.nstates;
if nargin == 4
    states2plot = varargin{1};
end

%%% quantity of interest
relstateerr = model.(whichindex).relstateerr;

%%% potential time unit conversion
unit = model.unit;
t_ref = unit.graphic.transftime(model.t_ref);

%%% graphical output
figure(figNr); clf
for k = 1:length(states2plot)
    plot(t_ref,relstateerr(:,states2plot(k)),'LineWidth',2,'Color',state2col(states2plot(k),:));
    hold on;
end
xlabel(sprintf('t [%s]',unit.graphic.time)); ylabel('rel state error'); title(sprintf('rel state error for %s-index',whichindex))
set(gca,'yscale','log'); legend(I.nmstatelegend(states2plot),'Location','eastoutside');
xlim(model.unit.graphic.xlim); 
model.makeplotbold(figNr); drawnow;

end

% -------------------------------------------------------------------------
function [] = plotnindices4state(figNr,I_state,model)

%%% potential unit conversion
unit = model.unit;
t_ref = unit.graphic.transftime(model.t_ref);

figure(figNr); clf;
legendtext = {};

for whichindex = {'env','pss','cneg','pneg'} 

    %%% check for existence of index
    if isfield(model,whichindex{:})
        if strcmp(whichindex{:},'cneg')
            legendtext = [legendtext {'cneg'}];
        elseif strcmp(whichindex{:},'pneg')
            legendtext = [legendtext {'pneg'}];
        else
            legendtext = [legendtext whichindex];
        end
    else
        continue
    end

    indx = model.(whichindex{:});

    if ismember(whichindex,{'env','cneg'})
        plot(t_ref,indx.nindex(:,I_state),'-','Linewidth',2);
        hold on;
    else
        plot(t_ref,indx.nindex(:,I_state),'--','Linewidth',2);
    end
end
% plot threshold for normalized index 
yline(model.threshold,'k--');

xlabel(sprintf('t [%s]',unit.graphic.time)); ylabel('normalised index');
title(sprintf('all indices for %s',model.I.nmstatelegend{I_state}))
%set(gca,'yscale','linear'); ylim([-0.05 1.05]); 
set(gca,'yscale','log'); ylim([1e-5 1e2]); 
xlim(model.unit.graphic.xlim);
legend([legendtext 'threshold'],'Location','eastoutside');
model.makeplotbold(figNr); drawnow;

end

% -------------------------------------------------------------------------
function [] = plotrelstateerr4state(figNr,I_state,model)


%%% potential unit conversion
unit = model.unit;
t_ref = unit.graphic.transftime(model.t_ref);

figure(figNr); clf;
legendtext = {};

for whichindex = {'env','pss','cneg','pneg'} 

    %%% check for existence of index
    if isfield(model,whichindex{:})
        if strcmp(whichindex{:},'cneg')
            legendtext = [legendtext {'cneg'}];
        elseif strcmp(whichindex{:},'pneg')
            legendtext = [legendtext {'pneg'}];
        else
            legendtext = [legendtext whichindex];
        end
    else
        continue
    end

    indx = model.(whichindex{:});

    if ismember(whichindex,{'env','cneg'})
        plot(t_ref,indx.relstateerr(:,I_state),'-','Linewidth',2);
        hold on;
    else
        plot(t_ref,indx.relstateerr(:,I_state),'--','Linewidth',2);
    end

end
% plot threshold for normalized index 
yline(model.threshold,'k--');

xlabel(sprintf('t [%s]',unit.graphic.time)); ylabel('rel state error');
title(sprintf('error for state %s',model.I.nmstatelegend{I_state}))
%set(gca,'yscale','linear'); ylim([-0.05 1.05]); 
set(gca,'yscale','log'); ylim([1e-5 1e2]); 
xlim(model.unit.graphic.xlim);
legend([legendtext 'threshold'],'Location','eastoutside');
model.makeplotbold(figNr); drawnow;

end

% -------------------------------------------------------------------------
function [] = plot_number_nir_above_threshold(figNr,model)

nnindexabovethreshold = sum([model.ir.nindex >= model.threshold],2);

%%% potential time unit conversion
unit = model.unit;
t_ref = unit.graphic.transftime(model.t_ref);

figure(figNr); clf
plot(t_ref,nnindexabovethreshold,'LineWidth',2);
xlabel(sprintf('t [%s]',unit.graphic.time)); ylabel(sprintf('number of nir indices above %d %%',100*model.threshold)); title(sprintf('sum(nir-index) >= %d/%',100*model.threshold))
set(gca,'yscale','lin'); 
% zoom in/out on the Y axis by changing the command below
xlim(model.unit.graphic.xlim); ylim([min(nnindexabovethreshold)-1 max(nnindexabovethreshold)+1])
model.makeplotbold(figNr); drawnow;
legend off

end


