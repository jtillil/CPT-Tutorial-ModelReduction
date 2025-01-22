function X0 = InitialCon(I)
%
%==========================================================================
%-- sources:
%  
% Andrea Y. Weiße, Diego A. Oyarzún, Vincent Danos, Peter S. Swain:
% 
% "Cellular trade-offs, gene expression, and growth"
% Proceedings of the National Academy of Sciences Mar 2015
% DOI: 10.1073/pnas.1416533112
%==========================================================================

% define initial state variable values
X0 = zeros(I.nstates,1);

X0(I.si)    = 0;         % intracellular nutrients
X0(I.a)     = 1000;        % a

X0(I.mm)    = 0;         % mRNAmetab
X0(I.mt)    = 0;         % mRNAt
X0(I.mq)    = 0;         % mRNAq
X0(I.mr)    = 0;         % mRNAr

X0(I.cm)   = 0;         % ribosomes + mRNAmetab
X0(I.ct)   = 0;         % ribosomes + mRNAt
X0(I.cq)   = 0;         % ribosomes + mRNAq 
X0(I.cr)   = 0;         % ribosomes + mRNAr

X0(I.zmm) = 0;       % strept + rib + mRNAmetab
X0(I.zmt) = 0;       % strept + rib + mRNAt
X0(I.zmq) = 0;       % strept + rib + mRNAq
X0(I.zmr) = 0;       % strept + rib + mRNAr

X0(I.em)    = 0;         % Metab. Enzyme
X0(I.et)    = 0;         % transporter enzyme
X0(I.q)     = 0;         % house queeping gene
X0(I.r)     = 10;        % free ribosomes 

end



