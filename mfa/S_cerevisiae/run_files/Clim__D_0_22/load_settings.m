%%% ADD DIR TO PATHS %%%
path_scripts = '../../../../mfa_scripts/';
path_runscripts = '../../../../run_scripts/';
path_case = './';
path_expmt = './';

gem_file = './scGEM_Clim__D_0_22_completeMFA.xlsx';
amm_file = './scAMM_Clim__D_0_22_completeMFA.xlsx';
expmt_files = {'data_expmt1.xlsx', 'data_expmt2.xlsx'};
expmt_ids = {'expmt1', 'expmt2'};

addpath(path_scripts, path_case, path_expmt, path_runscripts)

%%% SETTINGS %%%
% Scaling factor for random initial point (glucose uptake rate recommended) 
% and numbers of repeat initial start
randscale = 3.71;
repeat = 200;

% Constraining flux bounds for selected reactions
% vbs = {rxn_id, val, error}
% if vbfrombestfit = true, val will be replaced with
% the value of best fit extracted from res.mat
vbset = true;
vbfrombestfit = true;
vbs = {{'BIOMASS.f', 0, 0.0001},...
    {'EX_glc__D_e.f', 0, 0.0001}};

% Safe run on/off: if on, enabling saving with v7.3 setting
% Saved MATLAB objects are heavier but this is required for some
% run with large mapping model and/or labeling dataset
run_safe = false;

% Find best-fit enabling
run_bestfit = true;

% Run flux confidence interval estimation enabling
run_fconf = true;
