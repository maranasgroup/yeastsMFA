%% LOAD SETTINGS %%
run('load_settings.m');

%% START RUN %%
%% Find best fit %%
if run_bestfit
    if run_safe
        % Run safe enabled, save MATLAB objects with v7.3
        run_minSSR_func_safesave();
        
        % Partition heavy v7.3 MATLAB object for further processing
        load('res_v73.mat');
        save(strcat(path_expmt, 'res.mat'), 'model', 'mfamodel', 'emus',...
            'res', 'foptCell', 'residualCell');
        dAdv = emod.mdvsim.dAdv;
        save('emod.mdvsim.dAdv.mat', 'dAdv');
        dBdv = emod.mdvsim.dBdv;
        save('emod.mdvsim.dBdv.mat', 'dBdv');
        
    else
        run_minSSR_func();
    end
end

%% Find flux confidence interval %%
if run_fconf
    % Load best-fit result
    if run_safe
        load('res_v73.mat');
    else
        load('res.mat');
    end
    
    % Set bounds
    if vbset
        for j = 1:size(vbs,2)
            i_match = 0;
            vbmat = vbs{j};
            for i = 1:size(emod.vardata.flxdata,2)
                if strcmp(cell2mat(emod.vardata.flxdata(i).name), vbs{j}{1})
                    i_match = i;
                    break
                end
            end
            if vbfrombestfit
                vbs{j}{2} = res.fluxes(i_match).val;
            end
            emod.vardata.vb(i_match,1) = vbs{j}{2} - vbs{j}{3};
            emod.vardata.vb(i_match,2) = vbs{j}{2} + vbs{j}{3};
        end
    end
    
    % Determine flux confidence interval
    run_fluxconf_func(res,emod);
end