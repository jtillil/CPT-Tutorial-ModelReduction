function model = config2model(model, config)

%% check config
% check dimensions
if min(size(config)) ~= 1
    error(['Supplied config is not a vector. Size: ' num2str(size(config))])
end

% check length
if length(config) ~= model.I.nstates
    error(['nstates in config (' num2str(length(config)) ') does not equal nstates in I (' num2str(I.nstates) ')'])
end

% check that none of supplied configs are missing
if any(config == "")
    error(['The following states in config have not been assigned: ' num2str(find(config == ""))])
end

%% clean model config

model.I.dyn = [];
model.I.env = [];
model.I.neg = [];
model.I.pneg = [];
model.I.cneg = [];
model.I.pss = [];
model.I.qss = [];
model.I.tss = [];
model.I.irenv = [];
model.I.irenv_arith = [];
model.I.irenv_geom = [];
% model.I.con = [];

model.I.replaceODE = [];
model.I.replaceODEby = {};

%% write config to model

for k = 1:length(config)
    char_curr_conf = char(config(k));
    % fprintf(char_curr_conf)
    if char_curr_conf(1:3) ~= "con"
        model.I.(config(k)) = [model.I.(config(k)) k];
    else
        model.I.replaceODE = [model.I.replaceODE k];
        model.I.replaceODEby{end+1} = model.L.(model.L.nmconlaw{str2double(char_curr_conf(4:end))}).states;
    end
end

end

