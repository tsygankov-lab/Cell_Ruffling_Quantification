clear; clc;

root_fold = 'E:\Asef_Cdc42_Rac1_model\Rac1_regulator\2D_dynamic_cell_no_core\from_server_K1_edge_1.4_gamma_A_0.05_cell_R';

inp_folds = dir(root_fold);
inp_folds_ids = [inp_folds.isdir];
inp_folds = {inp_folds.name}';
inp_folds = inp_folds(inp_folds_ids);
inp_folds = inp_folds(3:end);

stp = 30;
M = 100;

for f_id = 1:length(inp_folds)
    inp_fold = inp_folds{f_id};
     
    info = imfinfo(fullfile(root_fold, inp_fold, 'trajectories', 'As.tif'));
    N_fr = length(info);
    
    MSK = cell(1,N_fr);
    
    for fr = 1:N_fr
        msk = double(imread(fullfile(root_fold, inp_fold, 'trajectories', 'As.tif'),fr)>0);
        msk = bwareaopen(msk, 10);
        msk = ~bwareaopen(~msk, 10);
        MSK{fr} = msk;
    end
    
    protr = zeros(N_fr-stp,1);
    retrc = zeros(N_fr-stp,1);
    
    for fr = 1:(N_fr-stp)
        m1 = MSK{fr};
        m2 = MSK{fr+stp};
        
        dif = m2 - m1;
        
        protr(fr) = sum(dif(:)>0);
        retrc(fr) = sum(dif(:)<0);
    end
    
    protr = (protr - mean(protr))/std(protr,1);
    retrc = (retrc - mean(retrc))/std(retrc,1);
    
    [protr_coeff, protr_fft, protr_rec] = fourier_shape_SH(protr,M);
    [retrc_coeff, retrc_fft, retrc_rec] = fourier_shape_SH(retrc,M);
    
    protr_y = protr_coeff;
    protr_y(protr_y==0) = NaN;
    
    retrc_y = retrc_coeff;
    retrc_y(retrc_y==0) = NaN;
    
    alpha_tr = 1;
    
    fig = figure('Position', [50 50 1800 900]);
    
    subplot(3,1,1);
    hold on;
    grid on;
    box on;
    plot(protr, 'Color', [0 1 0], 'LineWidth', 2, 'LineStyle', '-');
    plot(protr_rec, 'Color', [0 0 0], 'LineWidth', 1, 'LineStyle', '--');
    xlim([0, length(protr)]);
    title('Protrusion');
    legend({'origianl', 'approximation'}, 'Location', 'NorthWest');
    
    subplot(3,1,2);
    hold on;
    grid on;
    box on;
    plot(retrc, 'Color', [1 0 0], 'LineWidth', 2, 'LineStyle', '-');
    plot(retrc_rec, 'Color', [0 0 0], 'LineWidth', 1, 'LineStyle', '--');
    xlim([0, length(retrc)]);
    title('Retraction');
    legend({'origianl', 'approximation'}, 'Location', 'NorthWest');
    
    subplot(3,1,3);
    hold on;
    grid on;
    box on;
    scatter(1:M, protr_y, 'Marker', 'o', 'MarkerEdgeColor', [0 1 0], 'MarkerEdgeAlpha', alpha_tr, ...
        'MarkerFaceColor', [0 1 0], 'MarkerFaceAlpha', alpha_tr, 'SizeData', 50);
    scatter(1:M, retrc_y, 'Marker', 'o', 'MarkerEdgeColor', [1 0 0], 'MarkerEdgeAlpha', alpha_tr, ...
        'MarkerFaceColor', [1 0 0], 'MarkerFaceAlpha', alpha_tr, 'SizeData', 20);
    xticks(1:M);
    title('Fourier modes');
    legend({'protrusion','retraction'});
    saveas(fig, fullfile(root_fold, inp_fold, 'trajectories', strcat('fourier_modes_stp_', num2str(stp), '.png')));
    close(fig);
end
