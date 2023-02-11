clear; clc;

root_fold = 'E:\Asef_Cdc42_Rac1_model\ruffling_differentiator\2D_dynamic_cell';

folds = dir(root_fold);
ids = [folds.isdir];
folds = {folds(ids).name}';
folds = folds(3:end);

folds = {'cell_R_90_K1_edge_1.4_alpha_A_50_A_act_0.1_gamma_A_0.2'};

dp_vals = [2];

for f_id = 1:length(folds)
    fold = folds{f_id};
    disp(fold);
    fn = fullfile(root_fold, fold, 'trajectories', 'As.tif');
    for dp = dp_vals
        disp(dp);
        save_fold = fullfile(root_fold, fold, 'trajectories', strcat('dp_', num2str(dp)));
        mkdir(save_fold);
        
        try
            [elapsedTime,N_noINT,convg] = build_normals_YY2(fn,fn,dp,500,save_fold);
            elapsedTime = match_normals_YY2(fn,dp,save_fold);
            elapsedTime = quantify_along_normals_YY2(fn,fn,dp,save_fold);
        catch
            disp('problem with');
            disp(fn);
            disp(dp);
        end
    end
end
