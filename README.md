# CPT-Tutorial-ModelReduction

This repository contains supplemental code written in MATLAB to generate the results and figures of the CPT:PSP tutorial "Tackling High Dimensionality in QSP: Guiding Model Order Reduction with Index Analysis" and provides code and instructions to apply the presented approaches to your own QSP models.

Link to the publication: doi/...

## How to reproduce results and figures

To reproduce the results and figures presented in the tutorial, run the "MAIN.m" script located in the "/MAIN" folder while having the "/MAIN" folder as your working directory.

## How to apply the presented methods

Applying all steps of index-guided model order reduction to a large QSP model can take a considerable amount of time if run on an unsuitable machine. To test the capabilities, you can follow the following workflow with the small parallel pathways example model `modelSPP_with_crosstalk_full.mat` in the `/Core/modelfiles` folder.

### 1 Setup your QSP model

To be able to apply the presented model order reduction and index analysis framework, you first need to convert your QSP model into a compatible `model` MATLAB struct, requiring the following fields:
 - **model.name** [char]: name of the model
 - **model.scenario** [char]: name of the scenario in which to analyze the model
 - **model.I** [struct]: indexing struct describing the components (states, parameters) of the model, containing the following fields:
    - **model.I.nstates** [double]: scalar, number of states in the model
    - **model.I.nmstate** [cell]: 1-by-model.I.nstates cell array containing the state names as [char]
    - **model.I.<name_of_state>** [double]: for each state, a scalar indicating its location in model.X0, model.X_ref, and the input X to odefun and jacfun
    - **model.I.npar** [double]: scalar, number of parameters in the model
    - **model.I.nmpar** [cell]: 1-by-model.I.npar cell array containing the parameter names as [char]
    - **model.I.<name_of_parameter>** [double]: for each parameter, a scalar indicating its location in model.par and the input par to odefun and jacfun
 - **model.X0** [double]: model.I.nstates-by-1 vector containing the initial concentrations
 - **model.par** [double]: model.I.npar-by-1 vector containing the constant parameter values
 - **model.t_ref** [double]: ntsteps-by-1 vector indicating the timepoints at which the state trajectories should be evaluated
 - **model.X_ref** [double]: [also created by running create_minimal_model] ntsteps-by-model.I.nstates matrix containing the reference state trajectories
 - **model.odefun** [function_handle]: function with inputs (X [model.I.nstates-by-1 double containing the state concentrations at time *t*], par [model.I.npar-by-1 double containing the parameter values at time *t*]), returning dX [model.I.nstates-by-1 vector containing the derivative of the states at time *t*]
 - **model.jacfun** [function_handle]: function with inputs (X [model.I.nstates-by-1 double containing the state concentrations at time *t*], par [model.I.npar-by-1 double containing the parameter values at time *t*]), returning dF [model.I.nstates-by-model.I.nstates matrix containing the jacobian of the ODEs at time *t*], required for improved stability and efficiency of the ODE solvers

Then, run
```model = create_minimal_model(model)```
to automatically create all additionally required fields.

For a working example, load the saved model `modelBC_SV40_from_JKn_2024.mat` located in the `/Core/modelfiles` folder. It contains the blood coagulation network with the brown snake venom-fibrinogen relationship analyzed in this tutorial.

You will also be able to import your model from SBML (systems biology markup language) with the script `create_model_SBML.m` located in the `/Core/modelgeneration` folder. This capability is currently work in progress.

### 2 Apply index analysis

To apply index analysis to obtain insight into the model dynamics and guide the subsequent model order reduction, follow these steps:
 - Make the "/MAIN" folder your working directory.
 - Load your "model" struct into the workspace.
 - Run
   ```model = compute_and_analyse_indices_matlabfun(model)```
   to compute the ir and state-classification indices and relative state approximation errors to setup the `model` struct for further application of model order reduction. This step may take some time. It is recommended to perform this action on a modern multi-core CPU (or simply a server CPU, if available) and to set up the MATLAB parallel pool with a number of workers equivalent to half the number of logical cores available. On a 16-core, 32-thread AMD Ryzen 7950X CPU, computing all indices for the blood coagulation QSP model with 16 parallel workers took 6 hours and required 24GB of system memory.
 - The ir-indices are saved in the **model.ir** field alongside the controllability and observability indices in the **model.contr** and **model.obs** fields, the env, pneg, cneg, and pss state-classification indices are saved in the **model.env**, **model.pneg**, **model.cneg**, and **model.pss** fields, and the analysis of the indices is supplied in the **model.analysis** field.

### 3 Apply model order reduction

To apply index-guided model order reduction with the algorithm supplied in [Knöchel et. al. 2018], follow these steps:
 - Make the `/MAIN` folder your working directory.
 - Load your `model` struct into the workspace.
 - Run
   ```redmodel = mor_sequential_JKn_2018(model)```
   to reduce the model with the above-mentioned algorithm, using the sequence given by the ir-indices.
 - The reduction process and reduced model structure are documented in the console and saved in the newly created "redmodel" struct.

### 4 Compute ir-indices of the reduced model for validation

Validating if the reduced model retained the main mechanistic features of the original QSP model is a difficult task due to the extremely high dimensionality of the space of possible trajectories and can thus only been done locally. To do this via our presented index analysis approach, follow these steps:
 - Make the `/MAIN` folder your working directory.
 - Load your `redmodel` struct containing the reduced model obtained in step 3 into the workspace.
 - Run
   ```[ir, contr, obs] = compute_ir_indices_matlabfun(redmodel)```
   to compute the ir-indices for the reduced model.
 - Plot the (normalized) reduced indices against the original indices and compare if the dominant features are retained.
