function [L_old,newOutputIndex]=lumping_Aarons(model)
%%% lumping algorithm based on the paper by A. Dokoumetzidis and L.Aarons
%%% 'Proper lumping in systems biology models'
TOL             = model.lumping.TOL;
n               = length(model.Species);
P               = model.ODE.Matrix;
InitialValue    = model.InitialValue;
Time            = model.simTime;
Output          = model.Output;
[t,original_y]  = ode15s(@(t,C) model.ode(t,C,model),Time,InitialValue',model);
y_out           = original_y(:,Output);
L_old           = eye(n);
Error           = zeros(1,factorial(n)/(factorial(2)*factorial(n-2)));
newOutputIndex  = model.Output;

for i=1:n-1
    Combinations            = nchoosek(1:(n+1-i),2);
    L_pre                   = eye(n-i,n-i);
    L_poss                  = zeros(n-i,n-i+1,size(Combinations,1));
    L                       = zeros(n-i,n,size(Combinations,1));
    InvL                    = zeros(n,n-i,size(Combinations,1));
    model.ODE.SpecMatrix    = 'Lumping';
    
    for k=1:size(Combinations,1)
        L_poss(:,:,k)           = [L_pre(:,1:Combinations(k,2)-1),...
                                    L_pre(:,Combinations(k,1)),...
                                    L_pre(:,Combinations(k,2):end)];
        L(:,:,k)                = L_poss(:,:,k)*L_old;
        InvL(:,:,k)             = pinv(L(:,:,k));
        % if isfield(model,'environmentSpecies') 
        %     model.ODE.Matrix        = [L(:,:,k) zeros(size(L(:,:,k),1),length(model.environmentSpecies)); zeros(length(model.environmentSpecies),n) eye(length(model.environmentSpecies))]*P;
        %     model.ODE.InvLumpMatrix = P'*[InvL(:,:,k) zeros(n,length(model.environmentSpecies));zeros(length(model.environmentSpecies),size(L(:,:,k),1)) eye(length(model.environmentSpecies))];
        % else
        model.ODE.Matrix = L(:,:,k);
        model.ODE.InvLumpMatrix = InvL(:,:,k);
        % end
        if Combinations(k,:) < newOutputIndex
            newOutput(k) = newOutputIndex-1;
        else
            newOutput(k) = newOutputIndex;
        end
        switch model.lumping.flag
            case 1
                if Combinations(k,1)==newOutputIndex || Combinations(k,2)==newOutputIndex
                    Error(k)=Inf;
                continue;
                end
        end
        % if isfield(model,'environmentSpecies') 
        %     Initial = [L(:,:,k) zeros(size(L(:,:,k),1),length(model.environmentSpecies)); zeros(length(model.environmentSpecies),n) eye(length(model.environmentSpecies))]*InitialValue;
        % else
            Initial = L(:,:,k)*InitialValue';
        % end
        [~,y] = ode15s(@(t,C) model.OdeFcn(t,C,model),t,Initial,model);
        try
            Error(k) = relativeErrorL2(t,y_out,y(:,newOutput(k)));
        catch
           Error(k) = Inf;
        end
    end
    [~,Ind] = min(Error);
    if Combinations(Ind,:)<newOutputIndex
        newOutputIndex=newOutputIndex-1;
    end
    if min(Error)<TOL
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