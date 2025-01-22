function dX = ODEs(t,X,par,model)
%==========================================================================
%-- sources:
%  
% Andrea Y. Weiße, Diego A. Oyarzún, Vincent Danos, Peter S. Swain:
% 
% "Cellular trade-offs, gene expression, and growth"
% Proceedings of the National Academy of Sciences Mar 2015
% DOI: 10.1073/pnas.1416533112
%==========================================================================

I = model.I;

% initialise output vector
dX(size(X,1),1)= 0;

% assign variable names
a           = X(I.a);           % a
si          = X(I.si);          % interiorized nutrients

mm          = X(I.mm);          % mRNAm
mt          = X(I.mt);          % mRNAt
mq          = X(I.mq);          % mRNAq
mr          = X(I.mr);          % mRNAr

cm         = X(I.cm);         % ribosomes + mRNAmetab
ct         = X(I.ct);         % ribosomes + mRNAt
cq         = X(I.cq);         % ribosomes + mRNAq 
cr         = X(I.cr);         % ribosomes + mRNAr

zmm       = X(I.zmm);       % strept + rib + mRNAm
zmt       = X(I.zmt);       % strept + rib + mRNAt
zmq       = X(I.zmq);       % strept + rib + mRNAq
zmr       = X(I.zmr);       % strept + rib + mRNAr

em          = X(I.em);          % Metab. Enzyme
et          = X(I.et);          % transporter enzyme
q           = X(I.q);           % house queeping gene
r           = X(I.r);           % free ribosomes 


% assign parameter names
s0          = par(I.s0);       % External nutrients
ns          = par(I.ns);       % nutrient efficiency
vt          = par(I.vt);       % max nutrient import rate
vm          = par(I.vm);       % max enzymatic rate 
Km          = par(I.Km);       % enzymatic threshold
Kt          = par(I.Kt);       % Nutrient import threshold
thetar      = par(I.thetar);   % ribosome transcription threshold
thetax      = par(I.thetax);   % non-ribosomal transcription threshold
we          = par(I.we);       % max. enzyme transcription rate
wr          = par(I.wr);       % max. ribosome transcription rate
wq          = par(I.wq);       % max. q-transcription rate
Kq          = par(I.Kq);       % q-autoinhibition threshold
nq          = par(I.nq);       % q-autoinhibition Hill coeff. 
gmax        = par(I.gmax);     % max. transl. elongation rate 
Kgamma      = par(I.Kgamma);   % Net rate constant translational elongation (?)
nx          = par(I.nx);       % length of non-ribosomal proteins
nr          = par(I.nr);       % length of ribosomal proteins
kb          = par(I.kb);       % mRNA-ribosome binding rate
ku          = par(I.ku);       % mRNA-ribosome unbinding rate
dm          = par(I.dm);       % mRNA degradation rate 
Mref        = par(I.Mref);     % Cell reference mass 

cl          = par(I.cl);       % Chloramphenicol dose
k_cl        = par(I.k_cl);     % Chloramphenicol binding rate



% %%% import and metabolism rates
vimp        = et* vt*s0/(Kt+s0);
vcat        = em*vm*si/(Km+si);
 
% %%% transcription rates
Wm          = we * a/(thetax+a);
Wt          = we * a/(thetax+a);
Wq          = wq * a/(thetax+a)/(1+(q/Kq)^nq); 
Wr          = wr * a/(thetar+a);

% %%% translation rates
gamma      = gmax*a /(Kgamma+a); 
ttrate = (cq + cr + ct + cm)*gamma; 
 

 Vm         = gamma/nx * cm;
 Vt         = gamma/nx * ct;
 Vq         = gamma/nx * cq;
 Vr         = gamma/nr * cr;

% %%% dilution
 lam         = ttrate/Mref;

%------------------------- 

dX(I.si)= + vimp - vcat - lam*si; 
dX(I.a) = + ns*vcat - ttrate - lam*a;

dX(I.mm)= + Wm + ku*cm + Vm - kb*r*mm - dm*mm - lam*mm;
dX(I.mt)= + Wt + ku*ct + Vt - kb*r*mt - dm*mt - lam*mt;
dX(I.mq)= + Wq + ku*cq + Vq - kb*r*mq - dm*mq - lam*mq;
dX(I.mr)= + Wr + ku*cr + Vr - kb*r*mr - dm*mr - lam*mr;

dX(I.cm)= + kb*r*mm - ku*cm - Vm - k_cl*cl*cm - lam*cm; 
dX(I.ct)= + kb*r*mt - ku*ct - Vt - k_cl*cl*ct - lam*ct;
dX(I.cq)= + kb*r*mq - ku*cq - Vq - k_cl*cl*cq - lam*cq;
dX(I.cr)= + kb*r*mr - ku*cr - Vr - k_cl*cl*cr - lam*cr;

dX(I.zmm)= + k_cl*cl*cm - lam*zmm;
dX(I.zmt)= + k_cl*cl*ct - lam*zmt;
dX(I.zmq)= + k_cl*cl*cq - lam*zmq;
dX(I.zmr)= + k_cl*cl*cr - lam*zmr;

dX(I.em)= + Vm - lam*em;
dX(I.et)= + Vt - lam*et;
dX(I.q) = + Vq - lam*q;
dX(I.r) = + ku*cr + ku*ct + ku*cm + ku*cq +...
          Vr + Vr + Vt + Vm + Vq...
          - kb*r*mr - kb*r*mt - kb*r*mm - kb*r*mq -lam*r;

end 
