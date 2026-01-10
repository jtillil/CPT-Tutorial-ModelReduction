function config = I2config(I)

%% write config

config = repmat("", [1, I.nstates]);

config(I.dyn) = "dyn";
config(I.env) = "env";
config(I.pss) = "pss";
config(I.pneg) = "pneg";

end
