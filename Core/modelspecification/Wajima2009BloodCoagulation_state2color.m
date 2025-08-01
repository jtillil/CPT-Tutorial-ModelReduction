%%% Version: 19 Jun 2022
%%%
%%% call by: Wajima2009BloodCoagulation_state2color(I)
%%%
%%% This function assigns the colors to the given state variable to 
%%% facilitate plotting and to have always the same color for the specific 
%%% species.
%%%
%%% Citation:
%%% 
%%% Knoechel, Kloft and Huisinga, "Sensitivity based input-response index to 
%%% analyse and reduce large-scale signalling networks"
%%% PLOS Comp. Biology, 2020 (under review)
%%% 
%%% Authors: Jane Knoechel and Wilhelm Huisinga
%%%

function state2color = Wajima2009BloodCoagulation_state2color(I)

% define colors
color = [colormap('jet');colormap('HSV');colormap('Lines');colormap('jet');...
         colormap('HSV');colormap('Lines');colormap('jet');colormap('HSV');colormap('Lines')];
   
% define the state to color map     
map = struct('XII',color(1+(1-1)*7,:),'XIIa',color(2+(2-1)*7,:),'VIII',color(3+(3-1)*7,:),...
    'VIIIa',[0.301960796117783 0.745098054409027 0.933333337306976],'IX',color(5+(5-1)*7,:),'IXa',color(6+(6-1)*7,:),   ...
    'XI',color(7+(7-1)*7,:),'XIa',color(8+(8-1)*7,:),'VII',[0.6 0 0.2],...
    'VIIa',color(10+(10-1)*7,:),'X',color(11+(11-1)*7,:),'Xa',[0 0.498039215803146 0], ...
    'V', [0.678431391716003 0.921568632125854 1], 'Va',color(14+(14-1)*7,:),'Xa_Va',color(15+(15-1)*7,:),...
    'II',color(16+(16-1)*7,:),'IIa',color(17+(17-1)*7,:),'TAT',color(18+(18-1)*7,:), ...
    'Fg',color(19+(19-1)*7,:),'F',color(20+(20-1)*7,:),'XF',color(21+(21-1)*7,:),...
    'FDP',color(22+(22-1)*7,:),'D',color(23+(23-1)*7,:),'XIII',color(24+(24-1)*7,:), ...
    'XIIIa',color(25+(25-1)*7,:),'Pg',color(26+(26-1)*7,:),'P',color(27+(27-1)*7,:),...
    'PC',color(28+(28-1)*7,:),'APC',[0.466666668653488 0.674509823322296 0.18823529779911],'Tmod',color(30+(30-1)*7,:), ...
    'IIa_Tmod',color(31+(31-1)*7,:),'IXa_VIIIa',color(32+(32-1)*7,:),'TF',color(33+(33-1)*7,:),...
    'VII_TF',color(34+(34-1)*7,:),'VIIa_TF',[0.600000023841858 0.600000023841858 0],'TFPI',color(36+(36-1)*7,:), ...
    'Xa_TFPI',color(37+(37-1)*7,:),'VIIa_TF_Xa_TFPI',color(38+(38-1)*7,:),...
    'PS',color(39+(39-1)*7,:),'APC_PS',color(40+(40-1)*7,:),'Pk',color(41+(41-1)*7,:),...
    'K',color(42+(42-1)*7,:),'VK',color(43+(43-1)*7,:),'VKH2',color(44+(44-1)*7,:),...
    'VKO',color(45+(45-1)*7,:),'VK_p',color(46+(46-1)*7,:),'Awarf',color(47+(47-1)*7,:),...
    'Cwarf',color(48+(48-1)*7,:),'CA',color(49+(49-1)*7,:),'AEnox',color(50+(50-1)*7,:),...
    'AT_III_Heparin',color(51+(51-1)*7,:),'ENO_p',color(52+(52-1)*7,:),'AUC',color(53+(53-1)*7,:),...
    'AVenom',[1 0.843137254901961 0],'CVenom',color(55+(55-1)*7,:),'TaipanVenom',color(56+(56-1)*7,:),'ATIII',color(57+(57-1)*7,:),...
    'delayTaipan1',color(58+(58-1)*7,:),'delayTaipan2',color(59+(59-1)*7,:),'AVenom_Tiger',color(60+(60-1)*7,:),...
    'CVenom_Tiger',color(61+(61-1)*7,:),'AT_III_UFH',color(62+(62-1)*7,:));

% initialize output and assign colors
state2color = NaN(I.nstates,3);
for k = 1:I.nstates
    state2color(k,:) = map.(I.nmstate{k});
end

end