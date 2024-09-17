function config = model2config(model)

%% setup
I = model.I;
config = strings(1, I.nstates);

%% write to config
if isfield(I, 'dyn')
    if any(config(I.dyn) ~= "")
        error(['The following states in config have not been assigned: ' num2str(I.dyn(find(config(I.dyn) ~= "")))])
    end
    config(I.dyn) = "dyn";
end

if isfield(I, 'env')
    config(I.env) = "env";
end

if isfield(I, 'neg')
    config(I.neg) = "neg";
end

if isfield(I, 'pneg')
    config(I.pneg) = "pneg";
end

if isfield(I, 'cneg')
    config(I.cneg) = "cneg";
end

if isfield(I, 'qss')
    config(I.qss) = "qss";
end

if isfield(I, 'pss')
    config(I.pss) = "pss";
end

if isfield(I, 'irenv_arith')
    config(I.irenv_arith) = "irenv_arith";
end

if isfield(I, 'irenv_geom')
    config(I.irenv_geom) = "irenv_geom";
end

%% check config

% check length
if length(config) ~= I.nstates
    error(['nstates in config (' num2str(length(config)) ') does not equal nstates in I (' num2str(I.nstates) ')'])
end

% check that no configs from I missing
if any(config == "")
    error(['The following states in config have not been assigned: ' num2str(find(config == ""))])
end

end

