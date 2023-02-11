clear; clc;

root_fold = 'E:\Asef_Cdc42_Rac1_model\ruffling_differentiator\2D_dynamic_cell';

T = 10;
d = 10;

folds = dir(root_fold);

ids = [folds.isdir];
folds = {folds(ids).name}';
folds = folds(3:end);

folds = {'cell_R_90_K1_edge_1.4_alpha_A_50_A_act_0.1_gamma_A_0.2'};

for f_id = 1:length(folds)
    fold = folds{f_id};
    disp(fold);
    mkdir(fullfile(root_fold, fold, 'trajectories'));
    mkdir(fullfile(root_fold, fold, 'trajectories', 'mat'));
    
    N = length(dir(fullfile(root_fold, fold, 'mat', '*.mat')));
    
    for i = 1:T:N
        copyfile(fullfile(root_fold, fold, 'mat', strcat(num2str(i), '.mat')), ...
            fullfile(root_fold, fold, 'trajectories', 'mat', strcat(num2str(1+(i-1)/T), '.mat')));
    end
    
    files = dir(fullfile(root_fold, fold, 'trajectories', 'mat', '*.mat'));
    files = {files.name};
    N = length(files);
    
    V_vals = zeros(1, N);
    
    A_vals = zeros(1, N);
    C_vals = zeros(1, N);
    
    As_min_vals = zeros(1, N);
    As_max_vals = zeros(1, N);
    As_ampl_vals = zeros(1, N);
    
    A_min_vals = zeros(1, N);
    A_max_vals = zeros(1, N);
    A_ampl_vals = zeros(1, N);
    
    Cs_min_vals = zeros(1, N);
    Cs_max_vals = zeros(1, N);
    Cs_ampl_vals = zeros(1, N);
    
    C_min_vals = zeros(1, N);
    C_max_vals = zeros(1, N);
    C_ampl_vals = zeros(1, N);
    
    load(fullfile(root_fold, fold, 'trajectories', 'mat', strcat(num2str(1), '.mat')));
    S = size(Im_L);
    
    Im_L_all = zeros(S(1), S(2));
    
    for i = 1:N
        load(fullfile(root_fold, fold, 'trajectories', 'mat', strcat(num2str(i), '.mat')));
        Im = Im_L;
        As = As_L;
        A = A_L;
        Cs = Cs_L;
        C = C_L;
        V_vals(i) = sum(Im(:));
        A_vals(i) = sum(As(:)+A(:));
        C_vals(i) = sum(Cs(:)+C(:));
        As_min_vals(i) = min(As(Im==1));
        As_max_vals(i) = max(As(Im==1));
        As_ampl_vals(i) = max(As(Im==1))-min(As(Im==1));
        A_min_vals(i) = min(A(Im==1));
        A_max_vals(i) = max(A(Im==1));
        A_ampl_vals(i) = max(A(Im==1))-min(A(Im==1));
        Cs_min_vals(i) = min(Cs(Im==1));
        Cs_max_vals(i) = max(Cs(Im==1));
        Cs_ampl_vals(i) = max(Cs(Im==1))-min(Cs(Im==1));
        C_min_vals(i) = min(C(Im==1));
        C_max_vals(i) = max(C(Im==1));
        C_ampl_vals(i) = max(C(Im==1))-min(C(Im==1));
        
        Im_L_all = Im_L_all + Im_L;
    end
    
    ids_i = find(sum(Im_L_all,2));
    ids_j = find(sum(Im_L_all,1));
    
    i1 = ids_i(1);
    i2 = ids_i(end);
    j1 = ids_j(1);
    j2 = ids_j(end);
    
    save(fullfile(root_fold, fold, 'trajectories', 'data_inspection.mat'), ...
        'V_vals', 'A_vals', 'C_vals', 'As_min_vals', 'As_max_vals', ...
        'As_ampl_vals', 'A_min_vals', 'A_max_vals', 'A_ampl_vals', ...
        'Cs_min_vals', 'Cs_max_vals', 'Cs_ampl_vals', 'C_min_vals', ...
        'C_max_vals', 'C_ampl_vals', 'i1', 'i2', 'j1', 'j2');

    outp_file_As = fullfile(root_fold, fold, 'trajectories', 'As.tif');
    outp_file_Cs = fullfile(root_fold, fold, 'trajectories', 'Cs.tif');
    
    s = [i2-i1+2*d+1, j2-j1+2*d+1];
    A_data = zeros(s(1), s(2), N);
    As_data = zeros(s(1), s(2), N);
    C_data = zeros(s(1), s(2), N);
    Cs_data = zeros(s(1), s(2), N);
    
    A_min = min(A_min_vals);
    A_max = max(A_max_vals);
    As_min = min(As_min_vals);
    As_max = max(As_max_vals);
    C_min = min(C_min_vals);
    C_max = max(C_max_vals);
    Cs_min = min(Cs_min_vals);
    Cs_max = max(Cs_max_vals);
    
    for i = 1:N
        load(fullfile(root_fold, fold, 'trajectories', 'mat', strcat(num2str(i), '.mat')));
        A_data(:,:,i) = A_L(i1-d:i2+d, j1-d:j2+d);
        As_data(:,:,i) = As_L(i1-d:i2+d, j1-d:j2+d);
        C_data(:,:,i) = C_L(i1-d:i2+d, j1-d:j2+d);
        Cs_data(:,:,i) = Cs_L(i1-d:i2+d, j1-d:j2+d);
        
        As = As_L(i1-d:i2+d, j1-d:j2+d);
        As(As>0) = (As(As>0) - As_min)/(As_max-As_min);
        
        Cs = Cs_L(i1-d:i2+d, j1-d:j2+d);
        Cs(Cs>0) = (Cs(Cs>0) - Cs_min)/(Cs_max-Cs_min);
        if i == 1
            imwrite(As,outp_file_As);
            imwrite(Cs,outp_file_Cs);
        else
            imwrite(As,outp_file_As,'WriteMode','append');
            imwrite(Cs,outp_file_Cs,'WriteMode','append');
        end
    end
end

