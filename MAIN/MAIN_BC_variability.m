clear; clc;

% load model
load("modelBC_SV40_from_JKn_2024.mat")
config = repmat("dyn", [1, model.I.nstates]);
model.I = config2I(model.I, config, []);
par = model.par;
X0 = model.X0;

% generate virtual pop
seed = 1234;
rng(seed);
Npop = 100;

virtual_pop_par = zeros(Npop, length(par));
lognpdf_par = zeros(Npop, length(par));
for i = 1:length(par)
    if par(i) ~= 0
        m = par(i);
        v = (0.4*par(i))^2;
        mu = log(m^2 / sqrt(v + m^2));
        sigma = sqrt(log(v / (m^2) + 1));

        virtual_pop_par(:, i) = lognrnd(mu, sigma, [Npop 1]);
        lognpdf_par(:, i) = lognpdf(virtual_pop_par(:, i), mu, sigma);
    else
        lognpdf_par(:, i) = 1;
    end
end

virtual_pop_X0 = zeros(Npop, length(X0));
lognpdf_X0 = zeros(Npop, length(X0));
for i = 1:length(X0)
    if X0(i) ~= 0
        m = X0(i);
        v = (0.4*X0(i))^2;
        mu = log(m^2 / sqrt(v + m^2));
        sigma = sqrt(log(v / (m^2) + 1));

        virtual_pop_X0(:, i) = lognrnd(mu, sigma, [Npop 1]);
        lognpdf_X0(:, i) = lognpdf(virtual_pop_X0(:, i), mu, sigma);
    else
        lognpdf_X0(:, i) = 1;
    end
end
virtual_pop_X0(:, model.I.input) = X0(model.I.input);
lognpdf_X0(:, model.I.input) = 1;

% choose extreme individuals
lognpdf_par_indv = prod(lognpdf_par, 2);
[vals, idx] = sort(lognpdf_par_indv);
extreme_vals_par_5 = vals(1:5);
extreme_indv_par_5 = idx(1:5);

lognpdf_X0_indv = prod(lognpdf_X0, 2);
[vals, idx] = sort(lognpdf_X0_indv);
extreme_vals_X0_5 = vals(1:5);
extreme_indv_X0_5 = idx(1:5);

% compute indices for 5 extreme individuals
indv_nr = 2; % can be 1, 2, 3, 4, 5
indv_id = extreme_indv_X0_5(indv_nr);
model.X0 = virtual_pop_X0(indv_nr, :)';
[ir, contr, obs] = compute_ir_indices_matlabfun(model);
save(['results/irindices_BCSV40_from_JKn_full_extreme' char(string(indv_nr)) '_seed' char(string(seed)) '.mat'], "ir", "contr", "obs")