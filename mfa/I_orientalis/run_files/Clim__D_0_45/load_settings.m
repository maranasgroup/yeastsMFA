%%% ADD DIR TO PATHS %%%
path_scripts = '../../../../mfa_scripts/';
path_runscripts = '../../../../run_scripts/';
path_case = './';
path_expmt = './';

gem_file = './imIsor_GEM_completeMFA_Clim__D_0_45.xlsx';
amm_file = './imIsor_AMM_completeMFA_Clim__D_0_45.xlsx';
expmt_files = {'data_expmt1.xlsx', 'data_expmt2.xlsx'};
expmt_ids = {'expmt1', 'expmt2'};

addpath(path_scripts, path_case, path_expmt, path_runscripts)

%%% SETTINGS %%%
randscale = 6.18;
repeat = 100;
