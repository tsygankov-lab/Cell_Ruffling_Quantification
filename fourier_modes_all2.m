clear; clc;

root_fold = 'E:\Asef_Cdc42_Rac1_model\Rac1_regulator\2D_dynamic_cell_no_core\from_server_K1_edge_1.4_gamma_A_0.05_cell_R';

folds = dir(root_fold);
folds_ids = [folds.isdir];
folds = {folds.name}';
folds = folds(folds_ids);
folds = folds(3:end);

M = 100;
dp = 2;
t_step = 1;
R_disk = 2;

for f_id = 1:length(folds)
    fold = folds{f_id};
     
    load(fullfile(root_fold, fold, 'trajectories', strcat('dp_', num2str(dp)), strcat('tracked_dp', num2str(dp), '.mat')));
    N_frames = size(x_tr,1);
    prm = size(x_tr,2);
    file_location_AS = fullfile(root_fold, fold, 'trajectories', 'As.tif');
    AS_tr = read_concentration_kymograph(file_location_AS, N_frames, prm, t_step, R_disk, x_tr, y_tr);
    
    [AS_tr_coeff, AS_tr_fft, AS_tr_rec] = fourier_shape_SH2(AS_tr, M);
    [v_tr_coeff, v_tr_fft, v_tr_rec] = fourier_shape_SH2(v_tr, M);
    
    fig = figure('Position', [50 50 1800 900]);
    
    subplot(2,2,1);
    hold on;
    axis ij;
    imagesc(AS_tr);
    xlim([0, size(AS_tr,2)]);
    ylim([0, size(AS_tr,1)]);
    title('Original');
    
    subplot(2,2,2);
    hold on;
    axis ij;
    imagesc(AS_tr_rec);
    xlim([0, size(AS_tr_rec,2)]);
    ylim([0, size(AS_tr_rec,1)]);
    title('Approximation');
    
    subplot(2,2,[3,4]);
    hold on;
    grid on;
    box on;
    boxplot(AS_tr_coeff');
    ylim([0, max(AS_tr_coeff(:))]);
    title('Fourier coefficients');
    
    saveas(fig, fullfile(root_fold, fold, 'trajectories', strcat('dp_', num2str(dp)), 'AS_tr_fourier.png'));
    close(fig);
    
    fig = figure('Position', [50 50 1800 900]);
    
    subplot(2,2,1);
    hold on;
    axis ij;
    imagesc(v_tr);
    xlim([0, size(v_tr,2)]);
    ylim([0, size(v_tr,1)]);
    title('Original');
    
    subplot(2,2,2);
    hold on;
    axis ij;
    imagesc(v_tr_rec);
    xlim([0, size(v_tr_rec,2)]);
    ylim([0, size(v_tr_rec,1)]);
    title('Approximation');
    
    subplot(2,2,[3,4]);
    hold on;
    grid on;
    box on;
    boxplot(v_tr_coeff');
    ylim([0, max(v_tr_coeff(:))]);
    title('Fourier coefficients');
    
    saveas(fig, fullfile(root_fold, fold, 'trajectories', strcat('dp_', num2str(dp)), 'v_tr_fourier.png'));
    close(fig);
    
end
