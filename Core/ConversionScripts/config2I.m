function I = config2I(I, config, L)

%% check config
% check dimensions
if min(size(config)) ~= 1
    error(['Supplied config is not a vector. Size: ' num2str(size(config))])
end

% check length
if length(config) ~= I.nstates
    error(['nstates in config (' num2str(length(config)) ') does not equal nstates in I (' num2str(I.nstates) ')'])
end

% check that none of supplied configs are missing
if any(config == "")
    error(['The following states in config have not been assigned: ' num2str(find(config == ""))])
end

%% clean model config

I.dyn = [];
I.env = [];
I.neg = [];
I.pneg = [];
I.cneg = [];
I.pss = [];
I.qss = [];
I.tss = [];
I.irenv = [];
I.irenv_arith = [];
I.irenv_geom = [];
I.average = [];
I.mode = [];
I.ssenv = [];
I.constant = [];
I.constregr = [];
I.con = [];

I.replaceODE = [];
I.replaceODEby = {};

%% write config to model

for k = 1:length(config)
    char_curr_conf = char(config(k));
    % fprintf(char_curr_conf)
    if char_curr_conf(1:3) ~= "con"
        I.(config(k)) = [I.(config(k)) k];
    else
        if char_curr_conf(4) == "s"
            I.(config(k)) = [I.(config(k)) k];
        else
            I.con = [I.con k];
            I.replaceODE = [I.replaceODE k];
            I.replaceODEby{end+1} = L.(L.nmconlaw{str2double(char_curr_conf(4:end))}).states;
        end
    end
end

end
