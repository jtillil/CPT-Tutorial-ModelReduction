function config = I2config(I)

%% check I
% % check dimensions
% if min(size(config)) ~= 1
%     error(['Supplied config is not a vector. Size: ' num2str(size(config))])
% end
% 
% % check length
% if length(config) ~= I.nstates
%     error(['nstates in config (' num2str(length(config)) ') does not equal nstates in I (' num2str(I.nstates) ')'])
% end
% 
% % check that none of supplied configs are missing
% if any(config == "")
%     error(['The following states in config have not been assigned: ' num2str(find(config == ""))])
% end

%% write config

config = repmat("", [1, I.nstates]);

config(I.dyn) = "dyn";
config(I.env) = "env";
config(I.pss) = "pss";
config(I.pneg) = "pneg";
% I.cneg = [];
% config(I.dyn) = "dyn";
% I.pss = [];
% config(I.dyn) = "dyn";
% I.qss = [];
% config(I.dyn) = "dyn";
% I.tss = [];
% config(I.dyn) = "dyn";
% I.irenv = [];
% config(I.dyn) = "dyn";
% I.irenv_arith = [];
% config(I.dyn) = "dyn";
% I.irenv_geom = [];
% config(I.dyn) = "dyn";
% I.average = [];
% config(I.dyn) = "dyn";
% I.mode = [];
% config(I.dyn) = "dyn";
% I.ssenv = [];
% config(I.dyn) = "dyn";
% I.constant = [];
% config(I.dyn) = "dyn";
% I.constregr = [];
% config(I.dyn) = "dyn";
% I.con = [];
% config(I.dyn) = "dyn";
% 
% I.replaceODE = [];
% I.replaceODEby = {};


%% check resulting model config
% TODO
% has to account for states in conlaw

end
