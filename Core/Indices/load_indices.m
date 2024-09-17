%%% Version: 19 Jan 2023
%%%

function model = load_indices(model)

t_ref = model.t_ref; X_ref = model.X_ref;

no_consistent_data_found = 0;
for whichindex = model.I.nm_indices

    nmindx = whichindex{:};
    filename = [model.savenameroot '_' nmindx '_index.mat'];
    %%% check for existence of file name
    if ~exist(filename,'file')
        fprintf('\n No mat file for %s index found \n',whichindex{:});
        no_consistent_data_found = no_consistent_data_found +1;
        continue
    end    

    loadresults = load(filename);
    
    %%% check whether data is consistent with reference solution
    if ~( all(size(t_ref) == size(loadresults.model.t_ref)) && all(size(X_ref) == size(loadresults.model.X_ref)) && ...
          all(t_ref == loadresults.model.t_ref)  && all(X_ref(:) == loadresults.model.X_ref(:)) )
        fprintf('\n Mat file for %s index not consistent with reference model \n',nmindx);

        no_consistent_data_found = no_consistent_data_found +1;      
        continue
    end
        
    model.(nmindx) = loadresults.indx;
      
end

if length(model.I.nm_indices) == no_consistent_data_found
    fprintf('\n\n\n --> No consistent index data found at all!---PLEASE FIX! \n\n'); beep
    error(':-)');
end   

end


%%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%%% LOCAL SUB-ROUTINES
%%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

