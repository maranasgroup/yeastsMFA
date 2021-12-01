%%% ADD DIR TO PATHS %%%
path_scripts = '../../../../mfa_scripts/';
path_runscripts = '../../../../run_scripts/';
path_case = './';
path_expmt = './';

gem_file = './scGEM_completeMFA.xlsx';
amm_file = './scAMM_completeMFA.xlsx';
expmt_files = {'data_expmt1.xlsx', 'data_expmt2.xlsx'};
expmt_ids = {'expmt1', 'expmt2'};

addpath(path_scripts, path_case, path_expmt, path_runscripts)

%%% SETTINGS %%%
randscale = 25.55;
repeat = 100;
