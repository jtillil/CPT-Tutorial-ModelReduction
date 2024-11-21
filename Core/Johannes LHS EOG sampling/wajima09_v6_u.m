function [UI,DI,VI,idx,idp] = wajima09_v6_u(scenario,varargin)
%%% project: wajima09
%%% version: 0.6 (2022-02-01)
%%% authors: Undine Falkenhagen (0000-0003-0399-5886), Christian Himpe (0000-0003-2194-6754)
%%% license: BSD-2-Clause (opensource.org/licenses/BSD-2-Clause)
%%% summary: Blood coagulation model.

    global ODE;

    if (0 == nargin) || isempty(scenario), display('supply: ''invivo_snakevenom'', ''invivo_warfarin'', or ''invitro_pttest'''); return; end%if

    [idx,idp,xnames,pnames] = make_index();

    x0 = initial_state(idx);

    param = make_parameter(idp);

    switch scenario

        case 'invivo_snakevenom'

            dose_snake_venom = 0.0015; % [mg]
            SF_mg_to_nmol = 1e-3 / 2e5 * 1e9;
            u_ref = SF_mg_to_nmol * dose_snake_venom;
            x0(idx.AVenom) = x0(idx.AVenom) + u_ref;

            param(idp.v10) = 25000.0;
            param(idp.k10) = 1800.0;

            param(idp.v14) = 21000.0;
            param(idp.k14) = 30000.0;

            tspan = 0.0:1.0:40.0; % 40h

            C = sparse(1,idx.Fg,1.0,1,numel(fieldnames(idx)));
            g = @(x,u,p,t) C*x';

            invivo = 1.0;

        case 'invivo_warfarin'

            %u_ref = 4.0; % in [mg] unnecessary due to assumption of constant infusion instead of regular dose
            %x0(idx.Awarf) = x0(idx.Awarf) + u_ref;

            tspan = 0.0:24.0:10.0*24.0; % 30d

            %C = sparse(1,idx.VII,1.0,1,numel(fieldnames(idx)));
            
            C = sparse(1,idx.AUC,1.0,1,numel(fieldnames(idx)));
            g_in_vitro = @(x,u,p,t) t(find(C*x>5/12,1));
            u_in_vitro = sparse(1,idx.TF,100.0,1,numel(fieldnames(idx)));

            %g = @(x,u,p,t) ((x(idx.II,:).*x(idx.VII,:).*x(idx.X,:))./(x(idx.II,1)*x(idx.VII,1)*x(idx.X,1))).^-0.194;
            
            A_in_vitro = make_linear(idx,idp,param,x0,0.0);
            F_in_vitro = make_source(idx,idp,param,x0,0.0);
            make_Q_in_vitro = make_nonlinear(idx,idp,param,x0);
            xdot_in_vitro = @(t,x) A_in_vitro*x + F_in_vitro + make_Q_in_vitro(x)*x;

            g = @(x,u,p,t) arrayfun(@(n) response(0,g_in_vitro,[0.1/3600.0,30.0/3600.0],x(n,:)'./ 3.0,u_in_vitro,p),1:size(x,1));
            
            invivo = 1.0;

        case 'invitro_pttest'

            x0(idx.Tmod) = 0; % in the in vitro setting Tmod was assumed to be zero (Waj09)
            x0 = x0 ./ 3.0;   % and the dilution of the blood was considered by using 1/3 of the original initial values

            u_ref = 100.0; % in [nM]
            x0(idx.TF) = x0(idx.TF) + u_ref;

            tspan = 0.0:0.01/3600.0:30.0/3600.0; % 30s

            C = sparse(1,idx.AUC,1.0,1,numel(fieldnames(idx)));
            g = @(x,u,p,t) t(find(C*x'>5/12,1));

            invivo = 0.0;

        case 'invitro_pttest_low'

            x0(idx.Tmod) = 0; % in the in vitro setting Tmod was assumed to be zero (Waj09)
            x0 = x0 ./ 3.0;   % and the dilution of the blood was considered by using 1/3 of the original initial values

            u_ref = 5e-3; % in [nM]
            x0(idx.TF) = x0(idx.TF) + u_ref;

            tspan = 0.0:0.5/3600.0:240.0/3600.0; % 240s

            C = sparse(1,idx.AUC,1.0,1,numel(fieldnames(idx)));
            g = @(x,u,p,t) t(find(C*x'>5/12,1));

            invivo = 0.0;

    end%switch


    function [y,X] = wajima_solver(f,g,t,x0,u,p)
        persistent ctr;
        if isempty(ctr), ctr = 0; end%if

        A = make_linear(idx,idp,p,x0,invivo);

        F = make_source(idx,idp,p,x0,invivo);

        if strcmp(scenario,'invivo_warfarin')

            F(idx.Cwarf) = 0.017;
        end%if

        make_Q = make_nonlinear(idx,idp,p,x0);

        xdot = @(t,x) A*x + F + make_Q(x)*x;
        options.RelTol = 1e-5; options.AbsTol = 1e-7;
        [t,X] = ode15s(xdot,0:t(1):t(2),x0, options);

        y = g(X,0,p,t);

        fprintf('#');
        ctr = ctr + 1;
        if mod(ctr,40) == 0, fprintf('\n'); end%if

    end%function

    function [y,X] = response(f,g,t,x0,u,p)
        x0 = x0 + u';

        A_in_vitro = make_linear(idx,idp,p,x0,0.0);
        F_in_vitro = make_source(idx,idp,p,x0,0.0);
        make_Q_in_vitro = make_nonlinear(idx,idp,p,x0);
        xdot_in_vitro = @(t,x) A_in_vitro*x + F_in_vitro + make_Q_in_vitro(x)*x;

        [t,X] = ode15s(xdot_in_vitro,0:t(1):t(2),x0);

        y = g(X',0,p,t);

    end%function
    
    
    ODE = @wajima_solver;
    N = numel(x0);
%
    tic
    % with parameters
    W = emgr(@(x,u,p,t) 0,g,[0,N,1],tspan([2,end]),'i',param*[0.5,1.0,2],[1,0,0,0,1,0,0,0,3,0,1,0,0],[],[],x0,[],x0*[-0.5,1]);%,@(x,y) log(x + 1.0) * log(y + 1.0)); 
    %W = emgr(@(x,u,p,t) 0,g,[0,N,1],tspan([2,end]),'i',[param+log(0.5),param, param+log(2)],[1,0,0,0,1,0,0,0,3,2,1,0,0],[],[],x0,[],[log(0.5),log(2)]);    fprintf('\n');

    % only initial values
    %W = emgr(@(x,u,p,t) 0,g,[0,N,1],tspan([2,end]),'o',param,[1,0,0,0,1,0,0,0,3,0,0,0,0],[],[],x0,[],x0*[-0.5,1]); % 
    fprintf('\n');
    toc

    names = [xnames,pnames];

    %[UI,DI,VI] = svd(blkdiag(W{1},W{2}));
%     if scenario == "invivo_warfarin"
%         x_dyn = find(s_size<10^-8,1)-1;   % maybe actually take dyn(+qss) variables
%         WO = W(1:x_dyn, 1:x_dyn);       % Observability gramian
%         WM = W(1:x_dyn, x_dyn+1:end);   % Mixed block
%         WI = W(x_dyn+1:end, x_dyn+1:end) - (WM' * pinv(WO) * WM);
%         [UI,DI,VI] = svd(blkdiag(WO,WI));
%     else
        [UI,DI,VI] = svd(W);
%    end

    n=12;
    ui = abs(real(UI(:,1:n)));
    [s,o]=sort(ui*diag(DI(1:n,1:n)),'descend');% diag(DI(1:n,1:n))    ones(n,1)
    % with parameters
    figure; semilogy(1:numel(names),diag(DI)./max(diag(DI)),'linewidth',3); title('Observability');
    figure; h = bar3(ui(o(1:n+10),:)); yticks([1:numel(names)]); yticklabels(names(o(1:n+30)));
    % only initial values
    %figure; semilogy(1:numel(xnames),diag(DI)./max(diag(DI)),'linewidth',3); title('Observability');
    %figure; h = bar3(ui(o,:)); yticks([1:numel(xnames)]); yticklabels(xnames(o));
    
    %zcolor(h,ui);
    xlim([0,n+1]);
    ylim([0,n+30])%numel(names)]);
    

    

%{
     y = wajima_solver(@(x,u,p,t) 0,g,tspan([2,end]),x0,@(t) 0,param);
     figure; semilogy(y);
%}

end

function [idx,idp,xnames,pnames] = make_index()

    xnames = {'XII','XIIa','VIII','VIIIa','IX','IXa','XI','XIa','VII','VIIa', ...
              'X','Xa','V','Va','Xa_Va','II','IIa','TAT','Fg','F','XF','FDP', ...
              'D','XIII','XIIIa','Pg','P','PC','APC','Tmod','IIa_Tmod', ...
              'IXa_VIIIa','TF','VII_TF','VIIa_TF','TFPI','Xa_TFPI','VIIa_TF_Xa_TFPI', ...
              'PS','APC_PS','Pk','K','VK','VKH2','VKO','VK_p','Awarf','Cwarf', ...
              'CA','AEnox','AT_III_Heparin','ENO_p','AUC','AVenom','CVenom','TaipanVenom', ...
              'ATIII','delayTaipan1','delayTaipan2','AVenom_Tiger','CVenom_Tiger', ...
              'AT_III_UFH'};
          
    

    idx = struct();

    for k = 1:numel(xnames)
        idx = setfield(idx,xnames{k},k);
    end%for

    pnames = {'R1','R2',...
              'v1','v2','v3','v4','v5','v6','v7','v8','v9','v10','v11','v12','v13','v14','v15','v16','v17','v18','v19','v20',...
              'v21','v22','v23','v24','v25','v33','v34','v35','v36','v38','v39','v40','v41','v42','v43',...
              'k1','k2','k3','k4','k5','k6','k7','k8','k9','k10','k11','k12','k13','k14','k15','k16','k17','k18','k19','k20',...
              'k21','k22','k23','k24','k25','k33','k34','k35','k36','k38','k39','k40','k41','k42','k43',...
              'vtaipan','ktaipan',...
              'c26','c27','c28','c29','c30','c31','c32','c37','c44','c45',...
              ... %'lmax', ...
              'IC50',...
              'degXII','degXIIa','degVIII','degVIIIa','degIX','degIXa','degXI','degXIa','degVII','degVIIa','degX','degXa',...
              'degV','degVa','degXaVa','degII','degIIa','degTAT','degFg','degF','degXF','degFDP','degD','degXIII','degXIIIa',...
              'degPg','degP','degPC','degAPC','degTmod','degIIaTmod','degIXaVIIIa','degTF','degVIITF','degVIIaTF',...
              'degTFPI','degXaTFPI','degVIIaTFXaTFPI','degPS','degAPCPS','degPk','degK','degVK','degVK2','degCA',...
              'ka_Hep','Vc_Hep','Vp_Hep','Cl_Hep','Q_Hep','warf_dose','ka_Warf','Vd_Warf','Cl_Warf','VK_k12','VK_k21','VK_V',...
              'ka_Tiger','d_Tiger','ka_Brown','d_Brown','d_Taipan','ktrans_Taipan','inf_rate_UFH'};

    idp = struct();

    for k = 1:numel(pnames)
        idp = setfield(idp,pnames{k},k);
    end%for 
end

function x0 = initial_state(idx)

    N = numel(fieldnames(idx));
    x0 = zeros(N,1); 

    x0(idx.XII) = 375.0;
    x0(idx.XIIa) = 0.0;
    x0(idx.VIII) = 0.7;
    x0(idx.VIIIa) = 0.0;
    x0(idx.IX) = 89.6;
    x0(idx.IXa) = 0.0;
    x0(idx.XI) = 30.6;
    x0(idx.XIa) = 0.0;
    x0(idx.VII) = 10.0;
    x0(idx.VIIa) = 0.0;
    x0(idx.X) = 174.3;
    x0(idx.Xa) = 0.0;
    x0(idx.V) = 26.7;
    x0(idx.Va) = 0.0;
    x0(idx.Xa_Va) = 0.0;
    x0(idx.II) = 1394.4;
    x0(idx.IIa) = 0.0;
    x0(idx.TAT) = 0.0;
    x0(idx.Fg) = 8945.5;
    x0(idx.F) = 0.0;
    x0(idx.XF) = 0.0;
    x0(idx.FDP) = 0.0;
    x0(idx.D) = 0.0;
    x0(idx.XIII) = 70.3;
    x0(idx.XIIIa) = 0.0;
    x0(idx.Pg) = 2154.3;
    x0(idx.P) = 0.0;
    x0(idx.PC) = 60.0;
    x0(idx.APC) = 0.0;
    x0(idx.Tmod) = 50.0;
    x0(idx.IIa_Tmod) = 0.0;
    x0(idx.IXa_VIIIa) = 0.0;
    x0(idx.TF) = 0.0;
    x0(idx.VII_TF) = 0.0;
    x0(idx.VIIa_TF) = 0.0;
    x0(idx.TFPI) = 2.5;
    x0(idx.Xa_TFPI) = 0.0;
    x0(idx.VIIa_TF_Xa_TFPI) = 0.0;
    x0(idx.PS) = 300.0;
    x0(idx.APC_PS) = 0.0;
    x0(idx.Pk) = 450.0;
    x0(idx.K) = 0.0;
    x0(idx.VK) = 1.0;
    x0(idx.VKH2) = 0.1;
    x0(idx.VKO) = 0.1;
    x0(idx.VK_p) = 115.47541;
    x0(idx.Awarf) = 0.0;
    x0(idx.Cwarf) = 0.0;
    x0(idx.CA) = 0.0;
    x0(idx.AEnox) = 0.0;
    x0(idx.AT_III_Heparin) = 0.0; 
    x0(idx.ENO_p) = 0.0;
    x0(idx.AUC) = 0.0;
    x0(idx.AVenom) = 0.0;
    x0(idx.CVenom) = 0.0;
    x0(idx.TaipanVenom) = 0.0;
    x0(idx.ATIII) = 0.0;
    x0(idx.delayTaipan1) = 0.0;
    x0(idx.delayTaipan2) = 0.0;
    x0(idx.AVenom_Tiger) = 0.0;
    x0(idx.CVenom_Tiger) = 0.0;
    x0(idx.AT_III_UFH) = 0.0;
end

function p = make_parameter(idp)

    N = numel(fieldnames(idp));
    p = zeros(N,1);

    p(idp.R1) = 1.0 / 7.1;	% UFH
    p(idp.R2) = 1.0;		% UFH

    % V_max in [1/h]
    p(idp.v1)  = 50000.0;	% in VIII
    p(idp.v2)  = 50.0;		% in VIII
    p(idp.v3)  = 7.0;		% in IX
    p(idp.v4)  = 7.0;		% in XI
    p(idp.v5)  = 10.0;		% in XI
    p(idp.v6)  = 0.1;		% in VII
    p(idp.v7)  = 0.02;		% in X
    p(idp.v8)  = 2.0;		% in X
    p(idp.v9)  = 1.0;	% in X
    p(idp.v10) = 50000.0;	% in V
    p(idp.v11) = 50.0;		% in V
    p(idp.v12) = 100.0;	% in II
    p(idp.v13) = 9.0;		% in II
    p(idp.v14) = 20000.0;	% in I
    p(idp.v15) = 500.0;	% in I
    p(idp.v16) = 7.0;		% in I
    p(idp.v17) = 7.0;		% in I
    p(idp.v18) = 7.0;		% in I
    p(idp.v19) = 1.0;		% in I
    p(idp.v20) = 7.0;		% in XIII
    p(idp.v21) = 7.0;		% in Pg
    p(idp.v22) = 5.0;		% in Pg
    p(idp.v23) = 2.0;		% in Pg
    p(idp.v24) = 7.0;		% in PC
    p(idp.v25) = 2.0;		% in Xa_Va

    p(idp.v33) = 70.0;		% in TF
    p(idp.v34) = 900.0;	% in X
    p(idp.v35) = 70.0;		% in IX
    p(idp.v36) = 1000.0;	% in TF

    p(idp.v38) = 1.0;		% in VII
    p(idp.v39) = 1.0;		% in VII
    p(idp.v40) = 0.2;		% in VII
    p(idp.v41) = 7.0;		% in Hageman Factor XII
    p(idp.v42) = 70.0;		% in Hageman Factor XII
    p(idp.v43) = 0.0;		% in Pk | originally: 7.0, but: set either v43 or v42 to 0, otherwise unphysiologic activation of intrinsic pathway

    % K_m in [nM]
    p(idp.k1)  = 1.0;	% in VIII
    p(idp.k2)  = 1.0;		% in VIII
    p(idp.k3)  = 10.0;		% in IX
    p(idp.k4)  = 1.0;		% in XI
    p(idp.k5)  = 10.0;		% in XI
    p(idp.k6)  = 10.0;		% in VII
    p(idp.k7)  = 10.0;		% in X
    p(idp.k8)  = 0.1;		% in X
    p(idp.k9)  = 10.0;		% in X
    p(idp.k10) = 10.0;		% in V
    p(idp.k11) = 1.0;		% in V
    p(idp.k12) = 10.0;		% in II
    p(idp.k13) = 500.0;	% in II
    p(idp.k14) = 0.5;		% in I
    p(idp.k15) = 500.0;	% in I
    p(idp.k16) = 10.0;		% in I
    p(idp.k17) = 10.0;		% in I
    p(idp.k18) = 100.0;	% in I
    p(idp.k19) = 1.0;		% in I
    p(idp.k20) = 1.0;		% in XIII
    p(idp.k21) = 5000.0;	% in Pg
    p(idp.k22) = 10000.0;	% in Pg
    p(idp.k23) = 1.0;		% in Pg
    p(idp.k24) = 1.0;		% in PC
    p(idp.k25) = 1.0;		% in Xa_Va

    p(idp.k33) = 1.0;		% in TF
    p(idp.k34) = 200.0;	% in X
    p(idp.k35) = 1.0;		% in IX
    p(idp.k36) = 1.0;		% in TF

    p(idp.k38) = 10.0;		% in VII
    p(idp.k39) = 10.0;		% in VII
    p(idp.k40) = 10.0;		% in VII
    p(idp.k41) = 1.0;		% in Hageman Factor XII
    p(idp.k42) = 1.0;		% in Hageman Factor XII
    p(idp.k43) = 1.0;		% in Pk

    % taipan snake venom action on blood coagulation
    p(idp.vtaipan) = 70.0;
    p(idp.ktaipan) = 10.0;

    % c in [(nM*h)]
    p(idp.c26) = 0.01;		% in VIII and IXa_VIIIa
    p(idp.c27) = 0.5;		% in X
    p(idp.c28) = 0.5;		% in II and IIa_Tmod
    p(idp.c29) = 0.5;		% in VII
    p(idp.c30) = 0.1;		% in VII
    p(idp.c31) = 0.5;		% in TF and VIIa_TF_Xa_TFPI and Xa_TFPI
    p(idp.c32) = 0.5;		% in X and Xa_TFPI and TFPI

    p(idp.c37) = 0.5;		% in PC and PS and APC_PS

    p(idp.c44) = 0.85 * p(idp.R1);		% in II
    p(idp.c45) = 0.85;				% in X

    % warfarin PD parameters
    %p(idp.lmax) = 1.0;
    p(idp.IC50) = 0.34;

    % degradation in [1/h]
    p(idp.degXII) = 0.012;
    p(idp.degXIIa) = 20.0;
    p(idp.degVIII) = 0.058;
    p(idp.degVIIIa) = 20.0;
    p(idp.degIX) = 0.029;
    p(idp.degIXa) = 20.0;
    p(idp.degXI) = 0.10;
    p(idp.degXIa) = 20.0;
    p(idp.degVII) = 0.12;
    p(idp.degVIIa) = 20.0;
    p(idp.degX) = 0.018;
    p(idp.degXa) = 20.0;
    p(idp.degV) = 0.043;
    p(idp.degVa) = 20.0;
    p(idp.degXaVa) = 20.0;
    p(idp.degII) = 0.010;
    p(idp.degIIa) = 67.4;
    p(idp.degTAT) = 0.2;
    p(idp.degFg) = 0.032;
    p(idp.degF) = 0.050;
    p(idp.degXF) = 0.050;
    p(idp.degFDP) = 3.5;
    p(idp.degD) = 0.1;
    p(idp.degXIII) = 0.0036;
    p(idp.degXIIIa) = 0.69;
    p(idp.degPg) = 0.05;
    p(idp.degP) = 20.0;
    p(idp.degPC) = 0.050;
    p(idp.degAPC) = 20.4;
    p(idp.degTmod) = 0.050;
    p(idp.degIIaTmod) = 20.0;
    p(idp.degIXaVIIIa) = 20.0;
    p(idp.degTF) = 0.05;
    p(idp.degVIITF) = 0.7;
    p(idp.degVIIaTF) = 20.0;
    p(idp.degTFPI) = 20.0;
    p(idp.degXaTFPI) = 20.0;
    p(idp.degVIIaTFXaTFPI) = 20.0;
    p(idp.degPS) = 0.0165;
    p(idp.degAPCPS) = 20.0;
    p(idp.degPk) = 0.05;
    p(idp.degK) = 20.0;
    p(idp.degVK) = 0.2052;
    p(idp.degVK2) = 0.0228;
    p(idp.degCA) = 0.05;

    % heparin PK parameters (unfractioned heparin UFH or low-molecular weight heparin LMWH)
    p(idp.ka_Hep) = 0.255;
    p(idp.Vc_Hep) = 4.567;
    p(idp.Vp_Hep) = 29.6;
    p(idp.Cl_Hep) = 1.058;
    p(idp.Q_Hep) = 0.62;

    % warfarin PK parameters
    p(idp.warf_dose) = 4.0;
    p(idp.ka_Warf) = 1.0;
    p(idp.Vd_Warf) = 10.0;
    p(idp.Cl_Warf) = 0.2;

    % vitamin k parameters
    p(idp.VK_k12) = 0.0587;
    p(idp.VK_k21) = 0.0122;
    p(idp.VK_V) = 24.0;

    % tiger snake venom parameters
    p(idp.ka_Tiger) = 5.0;
    p(idp.d_Tiger) = 3.5;

    % brown snake venom parameters
    p(idp.ka_Brown) = 5.0;
    p(idp.d_Brown) = 3.5;

    % taipan snake venom parameters
    p(idp.d_Taipan) = 3.5;
    p(idp.ktrans_Taipan) = 0.99;

    p(idp.inf_rate_UFH) = 0;
end

function A = make_linear(idx,idp,p,x0,invivo)

    %
    degVKH2 = p(idp.degVK2) * x0(idx.VK) / x0(idx.VKH2);

    ke_Hep =  p(idp.Cl_Hep) / p(idp.Vc_Hep);
    k12_Hep = p(idp.Q_Hep) / p(idp.Vc_Hep);
    k21_Hep = p(idp.Q_Hep) / p(idp.Vp_Hep);

    ke_Warf = p(idp.Cl_Warf) / p(idp.Vd_Warf);

    N = numel(fieldnames(idx));
    A = sparse(N,N);

    % Diagonal Elements
    A(idx.XII,idx.XII) = -p(idp.degXII);
    A(idx.XIIa,idx.XIIa) = -p(idp.degXIIa);
    A(idx.VIII,idx.VIII) = -p(idp.degVIII);
    A(idx.VIIIa,idx.VIIIa) = -p(idp.degVIIIa);
    A(idx.IX,idx.IX) = -p(idp.degIX);
    A(idx.IXa,idx.IXa) = -p(idp.degIXa);
    A(idx.XI,idx.XI) = -p(idp.degXI);
    A(idx.XIa,idx.XIa) = -p(idp.degXIa);
    A(idx.VII,idx.VII) = -p(idp.degVII);
    A(idx.VIIa,idx.VIIa) = -p(idp.degVIIa);
    A(idx.X,idx.X) = -p(idp.degX);
    A(idx.Xa,idx.Xa) = -p(idp.degXa);
    A(idx.V,idx.V) = -p(idp.degV);
    A(idx.Va,idx.Va) = -p(idp.degVa);
    A(idx.Xa_Va,idx.Xa_Va) = -p(idp.degXaVa);
    A(idx.II,idx.II) = -p(idp.degII);
    A(idx.IIa,idx.IIa) = -p(idp.degIIa);
    A(idx.TAT,idx.TAT) = -p(idp.degTAT);
    A(idx.Fg,idx.Fg) = -p(idp.degFg);
    A(idx.F,idx.F) = -p(idp.degF);
    A(idx.XF,idx.XF) = -p(idp.degXF);
    A(idx.FDP,idx.FDP) = -p(idp.degFDP);
    A(idx.D,idx.D) = -p(idp.degD);
    A(idx.XIII,idx.XIII) = -p(idp.degXIII);
    A(idx.XIIIa,idx.XIIIa) = -p(idp.degXIIIa);
    A(idx.Pg,idx.Pg) = -p(idp.degPg);
    A(idx.P,idx.P) = -p(idp.degP);
    A(idx.PC,idx.PC) = -p(idp.degPC);
    A(idx.APC,idx.APC) = -p(idp.degAPC);
    A(idx.Tmod,idx.Tmod) = -p(idp.degTmod);
    A(idx.IIa_Tmod,idx.IIa_Tmod) = -p(idp.degIIaTmod);
    A(idx.IXa_VIIIa,idx.IXa_VIIIa) = -p(idp.degIXaVIIIa);
    A(idx.TF,idx.TF) = -p(idp.degTF);
    A(idx.VII_TF,idx.VII_TF) = -p(idp.degVIITF);
    A(idx.VIIa_TF,idx.VIIa_TF) = -p(idp.degVIIaTF);
    A(idx.TFPI,idx.TFPI) = -p(idp.degTFPI);
    A(idx.Xa_TFPI,idx.Xa_TFPI) = -p(idp.degXaTFPI);
    A(idx.VIIa_TF_Xa_TFPI,idx.VIIa_TF_Xa_TFPI) = -p(idp.degVIIaTFXaTFPI);
    A(idx.PS,idx.PS) = -p(idp.degPS);
    A(idx.APC_PS,idx.APC_PS) = -p(idp.degAPCPS);
    A(idx.Pk,idx.Pk) = -p(idp.degPk);
    A(idx.K,idx.K) = -p(idp.degK);
    A(idx.VK,idx.VK) = -p(idp.degVK) -p(idp.VK_k12);
    A(idx.VKH2,idx.VKH2) = -degVKH2;
    %A(idx.VKO,idx.VKO) = 0;
    A(idx.VK_p,idx.VK_p) = -p(idp.VK_k21);
    A(idx.Awarf,idx.Awarf) = -p(idp.ka_Warf);
    A(idx.Cwarf,idx.Cwarf) = -ke_Warf;
    A(idx.CA,idx.CA) = -p(idp.degCA);
    A(idx.AEnox,idx.AEnox) = -p(idp.ka_Hep);
    A(idx.AT_III_Heparin,idx.AT_III_Heparin) = -k12_Hep -ke_Hep;
    A(idx.ENO_p,idx.ENO_p) = -k21_Hep;
    %A(idx.AUC,idx.AUC) = 0;
    A(idx.AVenom,idx.AVenom) = -p(idp.ka_Brown);
    A(idx.CVenom,idx.CVenom) = -p(idp.d_Brown);
    A(idx.TaipanVenom,idx.TaipanVenom) = -p(idp.d_Taipan);
    A(idx.ATIII,idx.ATIII) = 0;
    A(idx.delayTaipan1,idx.delayTaipan1) = -p(idp.ktrans_Taipan);
    A(idx.delayTaipan2,idx.delayTaipan2) = -p(idp.ktrans_Taipan);
    A(idx.AVenom_Tiger,idx.AVenom_Tiger) = -p(idp.ka_Tiger);
    A(idx.CVenom_Tiger,idx.CVenom_Tiger) = -p(idp.d_Tiger);
    A(idx.AT_III_UFH,idx.AT_III_UFH) = ke_Hep; % positive entry is correct

    % Off-Diagonal Elements
    A(idx.IX,idx.VKH2) = p(idp.degIX) * x0(idx.IX) / x0(idx.VKH2) * invivo;
    A(idx.VII,idx.VKH2) = p(idp.degVII) * x0(idx.VII) / x0(idx.VKH2) * invivo;
    A(idx.X,idx.VKH2) = p(idp.degX) * x0(idx.X) / x0(idx.VKH2) * invivo;
    A(idx.II,idx.VKH2) = p(idp.degII) * x0(idx.II) / x0(idx.VKH2) * invivo;
    A(idx.PS,idx.VKH2) = p(idp.degPS) * x0(idx.PS) / x0(idx.VKH2) * invivo;
    A(idx.PC,idx.VKH2) = p(idp.degPC) * x0(idx.PC) / x0(idx.VKH2) * invivo;

    A(idx.AUC,idx.F) = 1.0;
    A(idx.TAT,idx.IIa) = p(idp.degIIa);
    A(idx.FDP,idx.Fg) = p(idp.degFg);
    A(idx.FDP,idx.F) = p(idp.degF);
    A(idx.D,idx.XF) = p(idp.degXF);
    A(idx.VK,idx.VK_p) = p(idp.VK_k21) / p(idp.VK_V);
    A(idx.VKO,idx.VKH2) = degVKH2;
    A(idx.VK_p,idx.VK) = p(idp.VK_k12) * p(idp.VK_V);
    A(idx.Cwarf,idx.Awarf) = p(idp.ka_Warf) / p(idp.Vd_Warf);
    A(idx.AT_III_Heparin,idx.AEnox) = p(idp.ka_Hep) / p(idp.Vc_Hep);
    A(idx.AT_III_Heparin,idx.ENO_p) = k21_Hep * p(idp.Vc_Hep);
    A(idx.ENO_p,idx.AT_III_Heparin) = k12_Hep * p(idp.Vc_Hep);
    A(idx.CVenom,idx.AVenom) = p(idp.ka_Brown);
    A(idx.CVenom_Tiger,idx.AVenom_Tiger) = p(idp.ka_Tiger);
    A(idx.delayTaipan1,idx.TaipanVenom) = p(idp.d_Taipan);
    A(idx.delayTaipan2,idx.delayTaipan1) = p(idp.ktrans_Taipan);
end

function F = make_source(idx,idp,p,x0,invivo)

    N = numel(fieldnames(idx));
    F = sparse(N,1);

    F(idx.XII) = p(idp.degXII) * x0(idx.XII) * invivo;
    F(idx.VIII) = p(idp.degVIII) * x0(idx.VIII) * invivo;
    F(idx.XI) = p(idp.degXI) * x0(idx.XI) * invivo;
    F(idx.V) = p(idp.degV) * x0(idx.V) * invivo;
    F(idx.Fg) = p(idp.degFg) * x0(idx.Fg) * invivo;
    F(idx.XIII) = p(idp.degXIII) * x0(idx.XIII) * invivo;
    F(idx.Pg) = p(idp.degPg) * x0(idx.Pg) * invivo;
    F(idx.Tmod) = p(idp.degTmod) * x0(idx.Tmod) * invivo;
    F(idx.TFPI) = p(idp.degTFPI) * x0(idx.TFPI) * invivo;
    F(idx.Pk) = p(idp.degPk) * x0(idx.Pk) * invivo;
    F(idx.VK) = p(idp.degVK) * x0(idx.VK) * invivo;

    F(idx.AT_III_UFH) = p(idp.inf_rate_UFH);
end

function y = f_vk(x,v,k)

    y = (v.*x) ./ (k+x);
end

function y = f_c(x,c)

    y = x ./ c;
end

function make_Q = make_nonlinear(idx,idp,p,x0)

    c46 = p(idp.c45) * p(idp.R2);

    degVKO = p(idp.degVK2) * x0(idx.VK) / x0(idx.VKO);

    N = numel(fieldnames(idx));

    % TODO preselect r, x and p components for bulk operation

    function Q = nonlin(x)

        r(1)  = f_vk(x(idx.IIa),p(idp.v1),0.000001*p(idp.k1));
        r(2)  = f_vk(x(idx.APC_PS),p(idp.v2),p(idp.k2));
        r(3)  = f_vk(x(idx.XIa),p(idp.v3),p(idp.k3));
        r(4)  = f_vk(x(idx.XIIa),p(idp.v4),p(idp.k4));
        r(5)  = f_vk(x(idx.IIa),p(idp.v5),p(idp.k5));
        r(6)  = f_vk(x(idx.IIa),p(idp.v6),p(idp.k6));
        r(7)  = f_vk(x(idx.IXa),p(idp.v7),p(idp.k7));
        r(8)  = f_vk(x(idx.IXa_VIIIa),p(idp.v8),p(idp.k8));
        r(9)  = f_vk(x(idx.VIIa),0.000000001*p(idp.v9),p(idp.k9));
        r(10) = f_vk(x(idx.IIa),p(idp.v10),p(idp.k10));
        r(11) = f_vk(x(idx.APC_PS),p(idp.v11),p(idp.k11));
        r(12) = f_vk(x(idx.Xa_Va)+x(idx.CVenom)+x(idx.TaipanVenom),p(idp.v12),p(idp.k12));
        r(13) = f_vk(x(idx.Xa)+x(idx.CVenom_Tiger),p(idp.v13),p(idp.k13));
        r(14) = f_vk(x(idx.IIa),p(idp.v14),p(idp.k14));
        r(15) = f_vk(x(idx.P),p(idp.v15),p(idp.k15));
        r(16) = f_vk(x(idx.XIIIa),p(idp.v16),p(idp.k16));
        r(17) = f_vk(x(idx.P),p(idp.v17),p(idp.k17));
        r(18) = f_vk(x(idx.P),p(idp.v18),p(idp.k18));
        r(19) = f_vk(x(idx.APC_PS),p(idp.v19),p(idp.k19));
        r(20) = f_vk(x(idx.IIa),p(idp.v20),p(idp.k20));
        r(21) = f_vk(x(idx.IIa),p(idp.v21),p(idp.k21));
        r(22) = f_vk(x(idx.F),p(idp.v22),p(idp.k22));
        r(23) = f_vk(x(idx.APC_PS),p(idp.v23),p(idp.k23));
        r(24) = f_vk(x(idx.IIa_Tmod),p(idp.v24),p(idp.k24));
        r(25) = f_vk(x(idx.APC_PS),p(idp.v25),p(idp.k25));

        r(26) = f_c(x(idx.IXa),p(idp.c26));
        r(27) = f_c(x(idx.Xa),p(idp.c27));
        r(28) = f_c(x(idx.IIa),p(idp.c28));
        r(29) = f_c(x(idx.TF),p(idp.c29));
        r(30) = f_c(x(idx.TF),p(idp.c30));
        r(31) = f_c(x(idx.VIIa_TF),p(idp.c31));
        r(32) = f_c(x(idx.Xa),p(idp.c32));

        r(33) = f_vk(x(idx.Xa),p(idp.v33),p(idp.k33));
        r(34) = f_vk(x(idx.VIIa_TF),p(idp.v34),p(idp.k34));
        r(35) = f_vk(x(idx.VIIa_TF),p(idp.v35),p(idp.k35));
        r(36) = f_vk(x(idx.TF),p(idp.v36),p(idp.k36));

        r(37) = f_c(x(idx.APC),p(idp.c37));

        r(38) = f_vk(x(idx.Xa),p(idp.v38),p(idp.k38));
        r(39) = f_vk(x(idx.VIIa_TF),p(idp.v39),p(idp.k39));
        r(40) = f_vk(x(idx.IXa),p(idp.v40),p(idp.k40));
        r(41) = f_vk(x(idx.CA),p(idp.v41),p(idp.k41));
        r(42) = f_vk(x(idx.K),p(idp.v42),p(idp.k42));
        r(43) = f_vk(x(idx.XIIa),p(idp.v43),p(idp.k43));

        r(44) = f_c(x(idx.IIa),p(idp.c44));
        r(45) = f_c(x(idx.Xa),p(idp.c45));
        r(46) = f_c(x(idx.IXa),c46);

        r(47) = (1.0 - f_vk(x(idx.Cwarf),1.0,p(idp.IC50))) * p(idp.degVK2); % p('lmax') -> 1.0
        r(48) = (1.0 - f_vk(x(idx.Cwarf),1.0,p(idp.IC50))) * degVKO;        % p('lmax') -> 1.0

        r(49) = f_vk(x(idx.delayTaipan2),p(idp.vtaipan),p(idp.ktaipan));

        Q = sparse(N,N);

        Q(idx.XII,idx.XII) = -r(41) -r(42);
        Q(idx.XIIa,idx.XII) = r(41) +r(42);
        Q(idx.VIII,idx.VIII) = -r(1);
        Q(idx.VIIIa,idx.VIII) = r(1);
        Q(idx.VIIIa,idx.VIIIa) = -r(2) -r(26);
        Q(idx.IX,idx.IX) = -r(35) -r(3);
        Q(idx.IXa,idx.IX) = r(35) +r(3);
        Q(idx.IXa,idx.AT_III_Heparin) = -r(46);
        Q(idx.IXa,idx.VIIIa) = -r(26);
        Q(idx.XI,idx.XI) = -r(5) -r(4);
        Q(idx.XIa,idx.XI) = r(5) +r(4);
        Q(idx.VII,idx.VII) = -r(6) -r(38) -r(39) -r(40) -r(49) -r(30);
        Q(idx.VIIa,idx.VII) = r(6) +r(38) +r(39) +r(40) +r(49);
        Q(idx.VIIa,idx.VIIa) = -r(29);
        Q(idx.X,idx.X) = -r(9) -r(34) -r(7) -r(8);
        Q(idx.Xa,idx.X) = r(9) +r(34) +r(7) +r(8);
        Q(idx.Xa,idx.AT_III_Heparin) = -r(45);
        Q(idx.Xa,idx.Va) = -r(27);
        Q(idx.Xa,idx.TFPI) = -r(32);
        Q(idx.V,idx.V) = -r(10);
        Q(idx.Va,idx.V) = r(10);
        Q(idx.Va,idx.Va) = -r(11) -r(27);
        Q(idx.Xa_Va,idx.Va) = r(27);
        Q(idx.Xa_Va,idx.Xa_Va) = -r(25);
        Q(idx.II,idx.II) = -r(12) -r(13);
        Q(idx.IIa,idx.II) = r(12) +r(13);
        Q(idx.IIa,idx.AT_III_Heparin) = -r(44);
        Q(idx.IIa,idx.Tmod) = -r(28);
        Q(idx.Fg,idx.Fg) = -r(14) -r(15);
        Q(idx.F,idx.Fg) = r(14);
        Q(idx.F,idx.F) = -r(16) -r(17);
        Q(idx.XF,idx.F) = r(16);
        Q(idx.XF,idx.XF) = -r(18) -r(19);
        Q(idx.FDP,idx.Fg) = r(15);
        Q(idx.FDP,idx.F) = r(17);
        Q(idx.D,idx.XF) = r(18) +r(19);
        Q(idx.XIII,idx.XIII) = -r(20);
        Q(idx.XIIIa,idx.XIII) = r(20);
        Q(idx.Pg,idx.Pg) = -r(21) -r(23) -r(22);
        Q(idx.P,idx.Pg) = r(21) +r(23) +r(22);
        Q(idx.PC,idx.PC) = -r(24);
        Q(idx.APC,idx.PC) = r(24);
        Q(idx.APC,idx.PS) = -r(37);
        Q(idx.VII_TF,idx.VII) = r(30);
        Q(idx.VII_TF,idx.VII_TF) = -r(36) -r(33);
        Q(idx.VII_TF,idx.Xa_TFPI) = -r(31);
        Q(idx.VIIa_TF,idx.VII_TF) = r(36) +r(33);
        Q(idx.VIIa_TF,idx.VIIa) = r(29);
        Q(idx.Pk,idx.Pk) = -r(43);
        Q(idx.K,idx.Pk) = r(43);
        Q(idx.VIIa_TF_Xa_TFPI,idx.Xa_TFPI) = r(31);
        Q(idx.Xa_TFPI,idx.TFPI) = r(32);
        Q(idx.Xa_TFPI,idx.Xa_TFPI) = -r(31);
        Q(idx.TFPI,idx.TFPI) = -r(32);
        Q(idx.TF,idx.VII) = -r(30);
        Q(idx.TF,idx.VIIa) = -r(29);
        Q(idx.Tmod,idx.Tmod) = -r(28);
        Q(idx.IIa_Tmod,idx.Tmod) = r(28);
        Q(idx.IXa_VIIIa,idx.VIIIa) = r(26);
        Q(idx.PS,idx.PS) = -r(37);
        Q(idx.APC_PS,idx.PS) = r(37);
        Q(idx.AT_III_Heparin,idx.AT_III_Heparin) = -r(44) -r(45) -r(46);
        Q(idx.AT_III_UFH,idx.AT_III_Heparin) = -r(44) -r(45) -r(46);
        Q(idx.VK,idx.VK) = -r(47);
        Q(idx.VK,idx.VKO) = r(48);
        Q(idx.VKH2,idx.VK) = r(47);
        Q(idx.VKO,idx.VKO) = -r(48);
    end%function

    make_Q = @nonlin;
end

function zcolor(handle,data)

    for k = 1:size(data,2)
        cdata = get(handle(k),'cdata');
        l = 1;
        for m = 0:6:(6*size(data,1)-6)
            cdata(m+1:m+6,:) = data(l,k);
            l = l+1;
        end%for
        set(handle(k),'cdata',cdata);
    end%for
end

