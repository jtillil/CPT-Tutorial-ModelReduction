function par = Par(I)
%==========================================================================
%-- sources:
%  
% Andrea Y. Weiße, Diego A. Oyarzún, Vincent Danos, Peter S. Swain:
% 
% "Cellular trade-offs, gene expression, and growth"
% Proceedings of the National Academy of Sciences Mar 2015
% DOI: 10.1073/pnas.1416533112
%==========================================================================

%%% Initializing parameter vector
par = NaN(I.npar,1);  
SF = 22000;

%%% external nutrient
par(I.s0)               = 1e4;

%%% nutrient import and metabolism
par(I.ns)               = 0.5;                  % Nutrient efficiency 
par(I.vt)               = 726.0;                % Max nutrient import rate [1/min]
par(I.vm)               = 5800.0;               % Max enzymatic rate [1/min]
par(I.Km)               = 1.0e3;                % Enzymatic threshold [molecs]
par(I.Kt)               = 1.0e3;                % Nutrient import threshold [molecs]

%%% transcription 
par(I.thetar)           = 426.8693338968694*SF;    % Ribosome transcription threshold [molecs]
par(I.thetax)           = 4.379733394834643*SF;    % Non-ribosomal transcription threshold [molecs]
par(I.we)               = 4.139172187824451;    % Max. enzyme transcription rate [molecs/min]
par(I.wr)               = 929.9678874564831;    % Max. ribosome transcription rate [molecs/min]
par(I.wq)               = 948.9349882947897;    % Max. q-transcription rate [molecs/min]

%%% autoinhibition for house keeping genes:
par(I.Kq)               = 1.522190403737490e+05;    % q-autoinhibition threshold [molecs]
par(I.nq)               = 4;                        % q-autoinhibition Hill coeff. 

%%% translation
par(I.gmax)             = 1260.0;               % Max. transl. elongation rate [aa/(min*molec)]
par(I.Kgamma)           = 7*SF;                    % Net rate constant translational elongation [molecs]
par(I.nx)               = 300.0;                % Length of non-ribosomal proteins [aa/molecs]
par(I.nr)               = 7549.0;               % Length of ribosomal proteins [aa/molecs]

par(I.kb)               = 1;                    % mRNA-ribosome binding rate [1/(min*molecs)]
par(I.ku)               = 1;                    % mRNA-ribosome unbinding rate [1/min]

par(I.dm)               = 0.1;                  % mRNA degradation rate [1/min]
par(I.Mref)             = 1.0e8;                % Cell reference mass [aa]

%%% chloramphenicol effect
par(I.cl)               = 0;                    % Chloramphenicol dose 
par(I.k_cl)             = 0.00599;              % Chloramphenicol binding rate [(min uM)^-1]

end
