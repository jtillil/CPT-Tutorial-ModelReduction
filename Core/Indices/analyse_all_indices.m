function analysis = analyse_all_indices(model)

fprintf('\n  This is analyse_all_indices ...\n')

threshold = model.threshold;

for whichindex = model.I.nm_indices

    nmindx = whichindex{:}; % name of index as a char string
    if ~isfield(model,nmindx)
        fprintf('\n NOTE: model structure does not contain %s index \n',nmindx)
        continue
    elseif strcmp(nmindx,'int') && ~isfield(model.int,'index') % field 'int' already exists
        continue
    end
    indx = model.(nmindx); 

    analysis.(nmindx) = [];
    switch nmindx
        case 'ir' %%% normalized ir-index
            
            max_nindex = max( indx.nindex(1:end-1,:),[],1 );
            [sorted_max_nindex,I_sorted_max_nindex] = sort( max_nindex,'descend' );

            %%% assign output
            %%%
            analysis.ir.max_nindex          = max_nindex;
            analysis.ir.sorted_max_nindex   = sorted_max_nindex;
            analysis.ir.I_sorted_max_nindex = I_sorted_max_nindex;

            I_sorted_max_nindex_above_threshold             = I_sorted_max_nindex( sorted_max_nindex >= threshold) ;
            analysis.ir.I_sorted_max_nindex_above_threshold = I_sorted_max_nindex_above_threshold;

            analysis.ir.plotsorted_max_nindex           = @(figNr) plotsorted_max_ir_nindex(figNr,analysis,model);
            analysis.ir.nmstates_above_nindex_threshold = model.I.nmstate( analysis.ir.I_sorted_max_nindex_above_threshold )';
            analysis.ir.nstates_above_nindex_threshold  = length( analysis.ir.nmstates_above_nindex_threshold );

            %%% plot sum of all indices
            analysis.ir.plotsum_ir_index = @(figNr) plotsum_ir_index(figNr,model);

        case {'env','pneg','cneg','pss'} %%% non-ir nindex and relative state errors

            %%% relative state errors
            %%% Note: last row is NaN by default, so exclude it from the analysis
            %%%
            max_relstateerr = max( indx.relstateerr(1:end-1,:),[],1,'includenan' );
            [sorted_max_relstateerr,I_max_relstateerr] = sort( max_relstateerr );
            I_sorted_max_relstateerr_below_threshold   = I_max_relstateerr( sorted_max_relstateerr < threshold) ;

            % assign output
            analysis.(nmindx).max_relstateerr        = max_relstateerr;
            analysis.(nmindx).sorted_max_relstateerr = sorted_max_relstateerr;
            analysis.(nmindx).I_sorted_max_relstateerr_below_threshold = I_sorted_max_relstateerr_below_threshold;

            %%% plot max of relative state error
            analysis.(nmindx).plotsorted_max_relstateerr = @(figNr,nmindx) plotsorted_max_relstateerr(figNr,analysis,nmindx,model);

            %%% non-ir nindex
            %%% Note: last row is NaN by default, so exclude it from the analysis)
            max_nindex = max( indx.nindex(1:end-1,:),[],1,'includenan' );
            [sorted_max_nindex,I_sorted_max_nindex] = sort( max_nindex );
            I_sorted_max_nindex_below_threshold     = I_sorted_max_nindex( sorted_max_nindex < threshold );

            %%% assign output
            %%%
            analysis.(nmindx).max_nindex          = max_nindex;
            analysis.(nmindx).sorted_max_nindex   = sorted_max_nindex;
            analysis.(nmindx).I_sorted_max_nindex = I_sorted_max_nindex;
            analysis.(nmindx).I_sorted_max_nindex_below_threshold = I_sorted_max_nindex_below_threshold;

            analysis.(nmindx).nmstates_below_nindex_threshold = model.I.nmstate( analysis.(nmindx).I_sorted_max_nindex_below_threshold )';
            analysis.(nmindx).nstates_below_nindex_threshold = length( analysis.(nmindx).nmstates_below_nindex_threshold ) ;

            analysis.(nmindx).I_states_below_nindex_and_state_threshold = intersect( I_sorted_max_nindex_below_threshold, I_sorted_max_relstateerr_below_threshold) ;
            analysis.(nmindx).nmstates_below_nindex_and_state_threshold = model.I.nmstate(analysis.(nmindx).I_states_below_nindex_and_state_threshold)';
            analysis.(nmindx).nstates_below_nindex_and_state_threshold  = length(analysis.(nmindx).nmstates_below_nindex_and_state_threshold);

            % assign plotting routines
            analysis.(nmindx).plotsorted_max_nindex = @(figNr) plotsorted_max_nindex(figNr,analysis,whichindex,model);

    end

end

end

%%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%%% LOCAL SUB-ROUTINES
%%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% -------------------------------------------------------------------------
function [] = plotsum_ir_index(figNr,model)

%%% potential unit conversion
unit = model.unit;
t_ref = unit.graphic.transftime(model.t_ref);

indx = model.ir.index;

figure(figNr); clf
plot(t_ref,sum(indx,2),'LineWidth',2);
title('sum of ir-index_k)')
xlabel(sprintf('t [%s]',unit.graphic.time));
legend('sum ir_k','Location','Northeast')
set(gca,'yscale','log');
xlim(model.setup.unit.graphic.xlim);

model.makeplotbold(figNr); drawnow;

end

% -------------------------------------------------------------------------
function [] = plotsorted_max_ir_nindex(figNr,analysis,model)

indx = analysis.ir;
threshold = model.threshold;

figure(figNr); clf
plot(indx.max_nindex(indx.I_sorted_max_nindex_above_threshold),'*','LineWidth',2,'Color',[0, 0.4470, 0.7410]);
hold on;
plot(indx.sorted_max_nindex,'*','LineWidth',2,'Color',[0.8500, 0.3250, 0.0980]);
yline(threshold,'b--');
plot(indx.max_nindex(indx.I_sorted_max_nindex_above_threshold),'*','LineWidth',2,'Color',[0, 0.4470, 0.7410]);
hold off;
title('sorted max(nir-index(k))')
xlabel('k'); ylabel('max(nir_k)')
legend('above threshold','below threshold',sprintf('threshold = %.2f',threshold),'Location','Northeast')
set(gca,'yscale','lin');
xlim(model.setup.unit.graphic.xlim);
model.makeplotbold(figNr); drawnow;

end

% -------------------------------------------------------------------------
function [] = plotsorted_max_nindex(figNr,analysis,nmindx,model)

indx = analysis.(nmindx);
threshold = model.threshold;

figure(figNr); clf
plot(indx.sorted_max_nindex,'*','LineWidth',2,'Color',[0.8500, 0.3250, 0.0980]);
hold on;
plot(indx.max_nindex(indx.I_sorted_max_nindex_below_threshold),'*','LineWidth',2,'Color',[0, 0.4470, 0.7410]);
yline(threshold,'b--');
hold off;
title(sprintf('sorted max(n%s-index(k))',nmindx))
legend(sprintf('max(n%s-index_k)',nmindx),sprintf('max(n%s-index_k)<threshold=%1.1e',nmindx,threshold),'Location','Southeast')
set(gca,'yscale','log');
xlim(model.setup.unit.graphic.xlim);
model.makeplotbold(figNr); drawnow;

end

% -------------------------------------------------------------------------
function [] = plotsorted_max_relstateerr(figNr,analysis,nmindx,model)

indx = analysis.(nmindx);
threshold = model.threshold;

figure(figNr); clf
plot(indx.max_relstateerr(indx.I_sorted_max_index),'*','LineWidth',2,'Color',[0.8500, 0.3250, 0.0980]);
hold on;
plot(indx.max_relstateerr(indx.I_sorted_max_index_below_threshold),'*','LineWidth',2,'Color',[0, 0.4470, 0.7410]);
yline(threshold,'b--');
hold off;
title(sprintf('max rel state err of sorted max %s index',nmindx))
legend(sprintf('max rel state err(%s-index_k)',nmindx),sprintf('max(%s-index_k)<threshold',nmindx),'Location','Southeast')
set(gca,'yscale','log');
xlim(model.setup.unit.graphic.xlim);
model.makeplotbold(figNr); drawnow;

end

