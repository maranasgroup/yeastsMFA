%%% ADD DIR TO PATHS %%%
path_scripts = '../../../../mfa_scripts/';
path_runscripts = '../../../../run_scripts/';
path_case = './';
path_expmt = './';

gem_file = './imIsor_GEM_SIMMER_completeMFA_Plim__D_0_34.xlsx';
amm_file = './imIsor_AMM_SIMMER_completeMFA_Plim.xlsx';
expmt_files = {'data_expmt1.xlsx', 'data_expmt2.xlsx'};
expmt_ids = {'expmt1', 'expmt2'};

addpath(path_scripts, path_case, path_expmt, path_runscripts)

%%% SETTINGS %%%
randscale = 4.06;
repeat = 100;
