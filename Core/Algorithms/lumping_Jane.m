function [L_old,newOutputIndex]=lumping_Jane(Model)
%%% lumping algorithm based on the paper by A. Dokoumetzidis and L.Aarons
%%% 'Proper lumping in systems biology models'
TOL             = Model.lumping.TOL;
n               = length(Model.Species);
P               = Model.ODE.Matrix;
InitialValue    = Model.InitialValue;
Time            = Model.simTime;
Output          = Model.Output;
[t,original_y]  = ode15s(@(t,C) Model.OdeFcn(t,C,Model),Time,InitialValue',Model);
y_out           = original_y(:,Output);
L_old           = eye(n);
Error           = zeros(1,factorial(n)/(factorial(2)*factorial(n-2)));
newOutputIndex  = Model.Output;

for i=1:n-1
    Combinations            = nchoosek(1:(n+1-i),2);
    L_pre                   = eye(n-i,n-i);
    L_poss                  = zeros(n-i,n-i+1,size(Combinations,1));
    L                       = zeros(n-i,n,size(Combinations,1));
    InvL                    = zeros(n,n-i,size(Combinations,1));
    Model.ODE.SpecMatrix    = 'Lumping';
    
    for k=1:size(Combinations,1)
        L_poss(:,:,k)           = [L_pre(:,1:Combinations(k,2)-1),...
                                    L_pre(:,Combinations(k,1)),...
                                    L_pre(:,Combinations(k,2):end)];
        L(:,:,k)                = L_poss(:,:,k)*L_old;
        InvL(:,:,k)             = pinv(L(:,:,k));
        if isfield(Model,'environmentSpecies') 
            Model.ODE.Matrix        = [L(:,:,k) zeros(size(L(:,:,k),1),length(Model.environmentSpecies)); zeros(length(Model.environmentSpecies),n) eye(length(Model.environmentSpecies))]*P;
            Model.ODE.InvLumpMatrix = P'*[InvL(:,:,k) zeros(n,length(Model.environmentSpecies));zeros(length(Model.environmentSpecies),size(L(:,:,k),1)) eye(length(Model.environmentSpecies))];
        else
            Model.ODE.Matrix = L(:,:,k);
            Model.ODE.InvLumpMatrix = InvL(:,:,k);
        end
        if Combinations(k,:)<newOutputIndex
            newOutput(k)    = newOutputIndex-1;
        else
            newOutput(k)   = newOutputIndex;
        end
        switch Model.lumping.flag
            case 1
                if Combinations(k,1)==newOutputIndex || Combinations(k,2)==newOutputIndex
                    Error(k)=Inf;
                continue;
                end
        end
        if isfield(Model,'environmentSpecies') 
            Initial = [L(:,:,k) zeros(size(L(:,:,k),1),length(Model.environmentSpecies)); zeros(length(Model.environmentSpecies),n) eye(length(Model.environmentSpecies))]*InitialValue;
        else
            Initial = L(:,:,k)*InitialValue';
        end
            [~,y]    = ode15s(@(t,C) Model.OdeFcn(t,C,Model),t,Initial,Model);
        try
            Error(k)        = relativeErrorL2(t,y_out,y(:,newOutput(k)));
        catch
           Error(k)        = Inf;
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