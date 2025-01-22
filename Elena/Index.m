function I = Index()
%==========================================================================
%-- sources:
%  
% Andrea Y. Weiße, Diego A. Oyarzún, Vincent Danos, Peter S. Swain:
% 
% "Cellular trade-offs, gene expression, and growth"
% Proceedings of the National Academy of Sciences Mar 2015
% DOI: 10.1073/pnas.1416533112
%==========================================================================

I.nmstate = {'si','a',...
             'mm','mt','mq','mr',...
             'cm', 'ct','cq', 'cr',... 
             'zmm', 'zmt', 'zmq', 'zmr',...
             'em', 'et', 'q', 'r'};


I.nstates  = length(I.nmstate);  % number (n) of states

for i = 1:I.nstates
    I.(I.nmstate{i}) = i;        % name to index mapping
end


%%% define parameter index structure
I.nmpar = {'s0','ns', 'vt', 'vm', 'Km', 'Kt', 'thetar',...
           'thetax', 'we', 'wr', 'wq', 'Kq', 'nq', 'gmax', 'Kgamma', 'nx',...
           'nr', 'kb', 'ku','dm','Mref',...
           'cl','k_cl'};


I.npar  = length(I.nmpar);       % number (n) of parameters

for i = 1:I.npar
    I.(I.nmpar{i}) = i;          % name to index mapping
end


%%% check for double naming and overlap between state and parameter names %%%

if I.nstates ~= length(unique(I.nmstate))
    fprintf('--> There are double named state names---PLEASE FIX!\n\n');
    error(':-)');
end
if I.npar ~= length(unique(I.nmpar))
    fprintf('--> There are double named parameter names---PLEASE FIX!\n\n');
    error(':-)');
end
if ~isempty(intersect(I.nmstate,I.nmpar))
    fprintf('--> There is an overlap betwenn state and parameter names :-( Please fix!\n\n');
    error('--> Ending here');
end


end