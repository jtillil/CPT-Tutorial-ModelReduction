
%%% ======================================================
% Economic Principles in Cell Biology (2024). 
% The Economic Cell Collective. No commercial publisher 
% Authors open access book. doi: 10.5281/zenodo.8156386.
% --------------------------------------------------------
% Chapter 8: Models of Growing cells
%%% ======================================================

% Problem 1.5 Simple coarse-grained models can generally be solved analytically.
% However, for models with higher level of granularity, like the one presented 
% in this section, reaching an analytical solution is highly complex. 
% Computational approaches that allow numerically solving high-dimensional
% systems are of great value.

%% ======================================================
% Exercise 1:  With the help of the provided code and following the detailed 
% description of the ODE system in Weisße et al. PNAS 2015, implement and 
% solve the system of ODEs. 
%
% The provided code is composed by 5 files that contain an incomplete 
% model implementation: 
% - Index: Maps each state variable and each parameter to a number for easy
% code readibility
% - Par: Creates an array that contains all parameter values stored by the 
% index structure 'I'.
% - InitialCon: Creates an array that contains the values of all initial
% conditions stored by the index structure 'I'.
% - ODEs: function where the ODE system is specified
% - ModelSim: This is the file used to simulate the model. It calls the 
% accompanying files and defines the options for the ode solver. 
% 
% Fill in the ODE system and solve it using the provided code:
clearvars;
% Initialize simulation----
% Import index  
  model.I = Index();
  I       = model.I;
% Parameter values
  model.par = Par(I);
% Initial values
  model.X0 = InitialCon(I);

  % model.tspan = [0, 1e6];
  model.tspan = [0, 4e3];

% Simulate model ---------
  model = ModelSim(model);
  t=model.t;
  X=model.X;

% Plot ---------------
  figure(1)
  plot(t,X,LineWidth=2 );
  ylabel('numbers of molecules')
  xlabel('Time [min]')
  legend(I.nmstate);



  figure(2)
  plot(model.t,model.X(:,I.a),LineWidth=2); hold on 
  plot(model.t,model.X(:,I.q),LineWidth=2);


%% REMOVE
model.X0 = model.X(end,:);
%% REMOVE
model.par(I.cl) = 0;
model = ModelSim(model);


plot(model.t,model.lam)


%% 
% Using this implementation, reproduce Monod s law, as seen in the inset of
% Figure 3

    s0_vec = [0:500:5e4]; % Define the vector of s0 to evaluate
    
    lam = []; %create a vector for storage
    for i = 1:length(s0_vec)
    model.par(I.s0) = s0_vec(i); % Modify s0
    model = ModelSim(model); % Simulate until steady state
    lam(i) = model.lam(end); % Save output (growth rate calculated in ModelSim.m)
    end 
    
    % Plot
    figure(2)
    plot(s0_vec,lam, LineWidth=2 )
    title("Monod's law")
    ylabel("growth rate [1/min]")
    xlabel("nutrient")


%% ========================================================================

% Exercise 2: The nutrient composition of the growth media is the main driver
% of increasing growth rates. Simulate the model to steady state for different 
% values of nutrient qualities. What model species are most impacted by an 
% increase in nutrient quality?

    ns_vec = [.1:.3:1]; % Define the vector of ns to evaluate
    
    ss_X = []; %create vector for storage
    for i = 1:length(ns_vec)
    model.par(I.ns) = ns_vec(i); % Modify ns parameter
    model = ModelSim(model);     % Simulate until steady state
    ss_X(:,i) = model.X(end,:)'; % Save output
    end 
    
    percChange = (ss_X./ss_X(:,1))*100; % Calculate percentage of change

    % Plot
    figure(3)
    plot(ns_vec,percChange, LineWidth=2)
    ylabel("Percentage of change")
    xlabel("nutrient quality")
    legend(I.nmstate);

    % Display ordered percentage of change
    totChange = table(percChange(:,4),'RowNames',I.nmstate, 'VariableNames',"% of change");
    disp(sortrows(totChange,1,'descend','MissingPlacement','last'));

% The species that are most affected by ns are: a>cr>r>...

%% ========================================================================

% Exercise 3: As seen in Figure 3, the addition of a drug that inhibits protein 
% synthesis results in an upregulation of the ribosomal fraction φR. Include 
% antibiotic function and reproduce Figure 3. How do the observed results relate 
% to your answer in question 2?


ns_vec = [0.08,0.11541599,0.16651064,0.24022489,0.3466,0.5]; % Define the vector of ns to evaluate
cl_vec = [0 2 4 8 12]; % Define the vector of chloramphenicol

% Vectors for storage 
Lam = [];
Rfr = [];
for i = 1:length(ns_vec)
    % Nutrient quality
    ns = ns_vec(i);
    model.par(I.ns) = ns;

    for j = 1:length(cl_vec)
        % Chlorampenicol dose
        cl = cl_vec(j);
        model.par(I.cl) = cl;
    
        % Simulate Model
         model = ModelSim(model);

         % Save output
         Lam(j,i) = model.lam(end);
         Rfr(j,i) = model.Rfr(end);
    end
end

% Plotting
figure(4)
plot(Lam, Rfr,'-x','linewidth',2); 
xlabel("Growth rate")
ylabel("Ribosome fraction")

