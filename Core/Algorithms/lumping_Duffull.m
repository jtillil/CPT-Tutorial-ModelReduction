function [lump_matrices, inv_lump_matrices, errors, out_states] = lumping_Duffull(model)
%%% lumping algorithm based on the paper by C. Hasegawa and S. B. Duffull
%%% 'Automated Scale Reduction of Nonlinear QSP Models'
%%% originally written by Jane Knöchel
%%% adapted by Johannes Tillil for use in the CPT Tutorial
%%% on index analysis and model reduction

%% setup
TOL             = 1e-1;
n               = model.I.nstates;
% P               = model.ODE.Matrix;
X0              = model.X0;
par             = model.par;
t_ref           = model.t_ref;
n_t             = length(t_ref);
output_state    = model.I.output;
[t_full,X_full] = ode15s(@(t,X) model.ode(t,X,par,model),t_ref,X0',model);
X_out           = X_full(:,output_state);
L_old           = eye(n);
Error           = zeros(1,factorial(n)/(factorial(2)*factorial(n-2)));
new_out_state   = model.I.output;

%% initial linearization
lin_TOL         = 1e-8;
lin_err         = 1;
f_y_t           = 1;
f               = zeros(1, n, n_t);
A_y_t           = 1;
A               = zeros(n, n, n_t);

% add initial f and A
for i = 1:n
    f(1, i, :) = X0(i);
    A(i, i, :) = X0(i);
end

% iteratively update f and A until tolerance met
while lin_err > lin_TOL
    
end

%% subsequent lumping
lump_matrices = cell(1, n-1);
inv_lump_matrices = cell(1, n-1);
errors = zeros(1, n-1);
out_states = zeros(1, n-1);

for i = 1:n-1
    Combinations            = nchoosek(1:(n+1-i),2);
    L_pre                   = eye(n-i,n-i);
    L_poss                  = zeros(n-i,n-i+1,size(Combinations,1));
    L                       = zeros(n-i,n,size(Combinations,1));
    InvL                    = zeros(n,n-i,size(Combinations,1));
    % model.ODE.SpecMatrix    = 'Lumping';
    
    for k = 1:size(Combinations,1)
        L_poss(:,:,k)           = [L_pre(:, 1:Combinations(k,2)-1), ...
                                   L_pre(:, Combinations(k,1)), ...
                                   L_pre(:, Combinations(k,2):end)];
        L(:,:,k)                = L_poss(:,:,k) * L_old;
        InvL(:,:,k)             = pinv(L(:,:,k));
        % if isfield(model,'environmentSpecies') 
        %     model.ODE.Matrix        = [L(:,:,k) zeros(size(L(:,:,k),1),length(model.environmentSpecies)); zeros(length(model.environmentSpecies),n) eye(length(model.environmentSpecies))]*P;
        %     model.ODE.InvLumpMatrix = P'*[InvL(:,:,k) zeros(n,length(model.environmentSpecies));zeros(length(model.environmentSpecies),size(L(:,:,k),1)) eye(length(model.environmentSpecies))];
        % else
        model.lumping.lumpmat = L(:,:,k);
        model.lumping.invlumpmat = InvL(:,:,k);
        % end
        if Combinations(k,:) < new_out_state
            new_out_states(k) = new_out_state-1;
        else
            new_out_states(k) = new_out_state;
        end
        % switch model.lumping.flag
        %     case 1
        %         if Combinations(k,1)==new_out_state || Combinations(k,2)==new_out_state
        %             Error(k)=Inf;
        %         continue;
        %         end
        % end
        % if isfield(model,'environmentSpecies') 
        %     X0_new = [L(:,:,k) zeros(size(L(:,:,k),1),length(model.environmentSpecies)); zeros(length(model.environmentSpecies),n) eye(length(model.environmentSpecies))]*X0;
        % else
            X0_new = L(:,:,k) * X0';
        % end
        [~,X_current] = ode15s(@(t,X) model.OdeFcn(t,X,model),t,X0_new,model);
        try
            Error(k) = relativeErrorL2(t_ref, X_full, X_current(:, new_out_states(k)));
        catch
           Error(k) = Inf;
        end
    end
    [~, Ind] = min(Error);
    if Combinations(Ind, :) < new_out_state
        new_out_state = new_out_state - 1;
    end
    if min(Error) < TOL
        L_old = L(:,:,Ind);
    else
        break;
    end
    clear Error
end
end

function Error=relativeErrorL2(t,X,Y)
Error=sqrt(trapz(t,(X-Y).^2)/trapz(t,X.^2));
end

function dX = ode_linearized(t, X, par, f, A)
    
end