clear; clc;

root_fold = 'E:\Asef_Cdc42_Rac1_model\ruffling_differentiator\2D_dynamic_cell';

folds = dir(root_fold);

ids = [folds.isdir];
folds = {folds(ids).name}';
folds = folds(3:end);

folds = {'cell_R_90_K1_edge_1.4_alpha_A_50_A_act_0.1_gamma_A_0.2'};

dp = 2;

w = 20;
t_step = 1;
R_disk = 2;

min_area = 10;

peak_mad = 2;
dip_mad = 2;

min_peaks_N = 5;
min_dips_N = 5;

c_N = 1000;
cm = jet(c_N);

% h = fspecial('gaussian', [10 1], 3);
% h = fspecial('gaussian', [12 1], 4);
% h = fspecial('gaussian', [15 1], 5);
% h = fspecial('gaussian', [15 1], 6);
% h = fspecial('gaussian', [5 1], 2);
% h = fspecial('gaussian', [5 1], 1);
h = fspecial('gaussian', [10 1], 4);
h = fspecial('gaussian', [20 1], 5);
h = fspecial('gaussian', [10 1], 2);

v_std_bin_N = 8;
v_step = 0.005;

v_nodes = -2:0.2:2;
v_w = 0.5;

rand_frac = 0.5;
v_dw = 0.05;
C_dw = 0.05;

std_A = 20;

diff_w = 1;

for fold_n = 1:length(folds)
    fold = folds{fold_n};
    disp(fold);
    mkdir(fullfile(root_fold, fold));
    tif_file_location = fullfile(root_fold, fold, 'trajectories', 'As.tif');
    save_fold = fullfile(root_fold, fold, 'trajectories', strcat('dp_', num2str(dp), '_UPD'));
    mkdir(save_fold);
    load(fullfile(root_fold, fold, 'trajectories', strcat('dp_', num2str(dp)), strcat('tracked_dp', num2str(dp), '.mat')));
    
    v_mean = mean(v_tr(:));
    v_std = std(v_tr(:));
    
    rm_trs = find(max(abs(v_tr),[],1)>std_A*v_std);
    sel_trs = find(max(abs(v_tr),[],1)<std_A*v_std);
    
    v_tr_sel = v_tr(:, sel_trs);
    x_tr_sel = x_tr(:, sel_trs);
    y_tr_sel = y_tr(:, sel_trs);
    
    fig = figure('Position', [50 50 1200 500]);
    
    subplot(1,2,1);
    hold on;
    axis image;
    axis ij;
    axis off;
    for p = 1:size(x_tr, 2)
        if ismember(p, rm_trs)
            plot(x_tr(:,p), y_tr(:,p), 'Color', [1 0 0], 'LineWidth', 0.5);
        else
            plot(x_tr(:,p), y_tr(:,p), 'Color', [0 0 0], 'LineWidth', 0.5);
        end
    end
    
    subplot(1,2,2);
    hold on;
    grid on;
    box on;
    histogram(v_tr_sel(:));
    plot([v_mean, v_mean], get(gca, 'YLim'), 'Color', [0 1 0], 'LineWidth', 1);
    plot([std_A*v_std, std_A*v_std], get(gca, 'YLim'), 'Color', [0 1 0], 'LineWidth', 1);
    plot([-std_A*v_std, -std_A*v_std], get(gca, 'YLim'), 'Color', [0 1 0], 'LineWidth', 1);
    saveas(fig, fullfile(save_fold, 'filtered_trajectories.png'));
    close(fig);
    
    v_tr = v_tr_sel;
    x_tr = x_tr_sel;
    y_tr = y_tr_sel;
    
    N_frames = size(x_tr,1);
    prm = size(x_tr,2);
    C_tr = read_concentration_kymograph(tif_file_location, N_frames, prm, t_step, R_disk, x_tr, y_tr);
    sC = imfilter(C_tr, h,'replicate');
    sVEL = imfilter(v_tr, h,'replicate');
    
    fig = figure('Position', [50 50 1800 500]);
    subplot(1,4,1);
    hold on;
    axis ij;
    colormap(jet);
    imagesc(C_tr);
    title('Biosensor Signal')
    xlabel('contour');
    ylabel('time');
    xlim([0, prm]);
    ylim([0, N_frames]);
    
    subplot(1,4,2);
    hold on;
    axis ij;
    colormap(jet);
    imagesc(sC);
    title('Smoothed Biosensor Signal')
    xlabel('contour');
    ylabel('time');
    xlim([0, prm]);
    ylim([0, N_frames]);
    
    subplot(1,4,3);
    hold on;
    axis ij;
    colormap(jet);
    imagesc(v_tr);
    title('Velocity')
    xlabel('contour');
    ylabel('time');
    xlim([0, prm]);
    ylim([0, N_frames]);
    
    subplot(1,4,4);
    hold on;
    axis ij;
    colormap(jet);
    imagesc(sVEL);
    title('Smoothed velocity')
    xlabel('contour');
    ylabel('time');
    xlim([0, prm]);
    ylim([0, N_frames]);
    saveas(fig, fullfile(save_fold, 'original_and_smoothed_kymographs.png'));
    close(fig);
    
    v_peak = median(sVEL(:)) + peak_mad*mad(sVEL(:));
    sVEL_peak_th = double(sVEL > v_peak);
    sVEL_peak_th = bwareaopen(sVEL_peak_th, min_area);
    peaks = sVEL_peak_th & islocalmax(sVEL,1);
    peaks(:,sum(peaks,1) < min_peaks_N) = 0;
    [I_peak, J_peak] = find(peaks);
    
    ids = (I_peak - w) > 0 & (I_peak + w) < size(v_tr,1);
    I_peak = I_peak(ids);
    J_peak = J_peak(ids);
    N_peak = length(I_peak);
    vw_peak_vals = zeros(N_peak, 2*w+1);
    Cw_peak_vals = zeros(N_peak, 2*w+1);
    for k = 1:N_peak
        vw_peak_vals(k,:) = v_tr(I_peak(k)-w:I_peak(k)+w, J_peak(k));
        Cw_peak_vals(k,:) = C_tr(I_peak(k)-w:I_peak(k)+w, J_peak(k));
    end
    xw_peak_vals = cumsum(vw_peak_vals,2);
    %dCw_peak_vals = diff(Cw_peak_vals,1,2);
    %dCw_peak_vals = Cw_peak_vals(:,1+diff_w:end) - Cw_peak_vals(:,1:end-diff_w);
    dCw_peak_vals = Cw_peak_vals(:,1+2*diff_w:end) - Cw_peak_vals(:,1:end-2*diff_w);
    
    v_dip = median(sVEL(:)) - dip_mad*mad(sVEL(:));
    sVEL_dip_th = double(sVEL < v_dip);
    sVEL_dip_th = bwareaopen(sVEL_dip_th, min_area);
    dips = sVEL_dip_th & islocalmin(sVEL,1);
    dips(:,sum(dips,1) < min_dips_N) = 0;
    [I_dip, J_dip] = find(dips);
    
    ids = (I_dip - w) > 0 & (I_dip + w) < size(v_tr,1);
    I_dip = I_dip(ids);
    J_dip = J_dip(ids);
    N_dip = length(I_dip);
    vw_dip_vals = zeros(N_dip, 2*w+1);
    Cw_dip_vals = zeros(N_dip, 2*w+1);
    for k = 1:N_dip
        vw_dip_vals(k,:) = v_tr(I_dip(k)-w:I_dip(k)+w, J_dip(k));
        Cw_dip_vals(k,:) = C_tr(I_dip(k)-w:I_dip(k)+w, J_dip(k));
    end
    xw_dip_vals = cumsum(vw_dip_vals,2);
    %dCw_dip_vals = diff(Cw_dip_vals,1,2);
    %dCw_dip_vals = Cw_dip_vals(:,1+diff_w:end) - Cw_dip_vals(:,1:end-diff_w);
    dCw_dip_vals = Cw_dip_vals(:,1+2*diff_w:end) - Cw_dip_vals(:,1:end-2*diff_w);
    
    vw_peak_SEM = std(vw_peak_vals)/sqrt(size(vw_peak_vals,1));
    vw_peak_ts = tinv([0.025  0.975], size(vw_peak_vals,1)-1);
    vw_peak_CI_1 = mean(vw_peak_vals,1) + vw_peak_ts(1)*vw_peak_SEM;
    vw_peak_CI_2 = mean(vw_peak_vals,1) + vw_peak_ts(2)*vw_peak_SEM;
    
    xw_peak_SEM = std(xw_peak_vals)/sqrt(size(xw_peak_vals,1));
    xw_peak_ts = tinv([0.025  0.975], size(xw_peak_vals,1)-1);
    xw_peak_CI_1 = mean(xw_peak_vals,1) + xw_peak_ts(1)*xw_peak_SEM;
    xw_peak_CI_2 = mean(xw_peak_vals,1) + xw_peak_ts(2)*xw_peak_SEM;
    
    Cw_peak_SEM = std(Cw_peak_vals)/sqrt(size(Cw_peak_vals,1));
    Cw_peak_ts = tinv([0.025  0.975], size(Cw_peak_vals,1)-1);
    Cw_peak_CI_1 = mean(Cw_peak_vals,1) + Cw_peak_ts(1)*Cw_peak_SEM;
    Cw_peak_CI_2 = mean(Cw_peak_vals,1) + Cw_peak_ts(2)*Cw_peak_SEM;
    
    dCw_peak_SEM = std(dCw_peak_vals)/sqrt(size(dCw_peak_vals,1));
    dCw_peak_ts = tinv([0.025  0.975], size(dCw_peak_vals,1)-1);
    dCw_peak_CI_1 = mean(dCw_peak_vals,1) + dCw_peak_ts(1)*dCw_peak_SEM;
    dCw_peak_CI_2 = mean(dCw_peak_vals,1) + dCw_peak_ts(2)*dCw_peak_SEM;
    
    vw_dip_SEM = std(vw_dip_vals)/sqrt(size(vw_dip_vals,1));
    vw_dip_ts = tinv([0.025  0.975], size(vw_dip_vals,1)-1);
    vw_dip_CI_1 = mean(vw_dip_vals,1) + vw_dip_ts(1)*vw_dip_SEM;
    vw_dip_CI_2 = mean(vw_dip_vals,1) + vw_dip_ts(2)*vw_dip_SEM;
    
    xw_dip_SEM = std(xw_dip_vals)/sqrt(size(xw_dip_vals,1));
    xw_dip_ts = tinv([0.025  0.975], size(xw_dip_vals,1)-1);
    xw_dip_CI_1 = mean(xw_dip_vals,1) + xw_dip_ts(1)*xw_dip_SEM;
    xw_dip_CI_2 = mean(xw_dip_vals,1) + xw_dip_ts(2)*xw_dip_SEM;
    
    Cw_dip_SEM = std(Cw_dip_vals)/sqrt(size(Cw_dip_vals,1));
    Cw_dip_ts = tinv([0.025  0.975], size(Cw_dip_vals,1)-1);
    Cw_dip_CI_1 = mean(Cw_dip_vals,1) + Cw_dip_ts(1)*Cw_dip_SEM;
    Cw_dip_CI_2 = mean(Cw_dip_vals,1) + Cw_dip_ts(2)*Cw_dip_SEM;
    
    dCw_dip_SEM = std(dCw_dip_vals)/sqrt(size(dCw_dip_vals,1));
    dCw_dip_ts = tinv([0.025  0.975], size(dCw_dip_vals,1)-1);
    dCw_dip_CI_1 = mean(dCw_dip_vals,1) + dCw_dip_ts(1)*dCw_dip_SEM;
    dCw_dip_CI_2 = mean(dCw_dip_vals,1) + dCw_dip_ts(2)*dCw_dip_SEM;
    
    corr_vel_C = zeros(2*w+1,prm);
    
    for p = 1:prm
        
        v_p = v_tr(w+1:N_frames-t_step-w,p);
        
        for t = 1:2*w+1
            C_pt = C_tr(t:N_frames-t_step-2*w-1+t,p);
            corr_vel_C(t,p) = corr(v_p, C_pt);
        end
    end
    
    [v_vals, ids] = sort(v_tr(:));
    C_vals = C_tr(:);
    C_vals = C_vals(ids);
    
    v_vals2 = -4:v_step:4;
    C_vals2 = zeros(size(v_vals2));
    C_stds2 = zeros(size(v_vals2));
    
    for i = 1:length(v_vals2)
        v = v_vals2(i);
        ids = (v_vals > v - v_step) & (v_vals < v + v_step);
        C_vals2(i) = mean(C_vals(ids));
        C_stds2(i) = std(C_vals(ids));
    end
    
    v_std = std(v_tr);
    v_std_min = min(v_std);
    v_std_max = max(v_std);
    
    v_tr_all_vals = v_tr(:);
    C_tr_all_vals = C_tr;
    ids = (v_tr_all_vals >= min(v_nodes)) & (v_tr_all_vals <= max(v_nodes));
    v_tr_all_vals = v_tr_all_vals(ids);
    C_tr_all_vals = C_tr_all_vals(ids);
    C_tr_all_vals = (C_tr_all_vals - mean(C_tr_all_vals))/std(C_tr_all_vals);
    C_tr_all_mean_vals = zeros(size(v_nodes));
    C_tr_all_CI_1_vals = zeros(size(v_nodes));
    C_tr_all_CI_2_vals = zeros(size(v_nodes));
    
    v_tr_peak_vals = vw_peak_vals(:);
    C_tr_peak_vals = Cw_peak_vals(:);
    ids = (v_tr_peak_vals >= min(v_nodes)) & (v_tr_peak_vals <= max(v_nodes));
    v_tr_peak_vals = v_tr_peak_vals(ids);
    C_tr_peak_vals = C_tr_peak_vals(ids);
    C_tr_peak_vals = (C_tr_peak_vals - mean(C_tr_peak_vals))/std(C_tr_peak_vals);
    C_tr_peak_mean_vals = zeros(size(v_nodes));
    C_tr_peak_CI_1_vals = zeros(size(v_nodes));
    C_tr_peak_CI_2_vals = zeros(size(v_nodes));
    
    v_tr_dip_vals = vw_dip_vals(:);
    C_tr_dip_vals = Cw_dip_vals(:);
    ids = (v_tr_dip_vals >= min(v_nodes)) & (v_tr_dip_vals <= max(v_nodes));
    v_tr_dip_vals = v_tr_dip_vals(ids);
    C_tr_dip_vals = C_tr_dip_vals(ids);
    C_tr_dip_vals = (C_tr_dip_vals - mean(C_tr_dip_vals))/std(C_tr_dip_vals);
    C_tr_dip_mean_vals = zeros(size(v_nodes));
    C_tr_dip_CI_1_vals = zeros(size(v_nodes));
    C_tr_dip_CI_2_vals = zeros(size(v_nodes));
    
    for i = 1:length(v_nodes)
        v = v_nodes(i);
        
        ids = (v_tr_all_vals >= v-v_w) & (v_tr_all_vals <= v+v_w);
        N_i = sum(ids);
        ts = tinv([0.025  0.975], N_i-1);
        
        C_i_mean = mean(C_tr_all_vals(ids));
        C_i_std = std(C_tr_all_vals(ids));
        C_tr_all_mean_vals(i) = C_i_mean;
        C_tr_all_CI_1_vals(i) = C_i_mean + ts(1)*C_i_std/sqrt(N_i);
        C_tr_all_CI_2_vals(i) = C_i_mean + ts(2)*C_i_std/sqrt(N_i);
        
        ids = (v_tr_peak_vals >= v-v_w) & (v_tr_peak_vals <= v+v_w);
        N_i = sum(ids);
        ts = tinv([0.025  0.975], N_i-1);
        
        C_i_mean = mean(C_tr_peak_vals(ids));
        C_i_std = std(C_tr_peak_vals(ids));
        C_tr_peak_mean_vals(i) = C_i_mean;
        C_tr_peak_CI_1_vals(i) = C_i_mean + ts(1)*C_i_std/sqrt(N_i);
        C_tr_peak_CI_2_vals(i) = C_i_mean + ts(2)*C_i_std/sqrt(N_i);
        
        ids = (v_tr_dip_vals >= v-v_w) & (v_tr_dip_vals <= v+v_w);
        N_i = sum(ids);
        ts = tinv([0.025  0.975], N_i-1);
        
        C_i_mean = mean(C_tr_dip_vals(ids));
        C_i_std = std(C_tr_dip_vals(ids));
        C_tr_dip_mean_vals(i) = C_i_mean;
        C_tr_dip_CI_1_vals(i) = C_i_mean + ts(1)*C_i_std/sqrt(N_i);
        C_tr_dip_CI_2_vals(i) = C_i_mean + ts(2)*C_i_std/sqrt(N_i);
    end
    
    fig = figure('Position', [50 50 800 800]);
    hold on;
    axis image; axis ij; axis off;
    for p = 1:size(x_tr, 2)
        c_i = fix((v_std(p)-v_std_min)/(v_std_max-v_std_min)*(c_N-1))+1;
        plot(x_tr(:,p), y_tr(:,p), 'Color', cm(c_i,:), 'LineWidth', 0.5);
    end
    title('v std for trajectories');
    saveas(fig, fullfile(save_fold, 'trajectories_velocity_std.png'));
    close(fig);
    
    fig = figure('Position', [50 50 1800 500]);
    subplot(1,4,1);
    hold on;
    axis ij;
    imagesc(sVEL);
    title('Velocity (smoothed)')
    xlabel('contour');
    ylabel('time');
    xlim([0, prm]);
    ylim([0, N_frames]);
    
    subplot(1,4,2);
    hold on;
    axis ij;
    imagesc(sVEL_peak_th);
    title('Velocity threshold')
    xlabel('contour');
    ylabel('time');
    xlim([0, prm]);
    ylim([0, N_frames]);
    
    subplot(1,4,3);
    hold on;
    axis ij;
    imagesc(v_tr);
    scatter(J_peak, I_peak, ...
        'Marker', 'o', 'SizeData', 5, ...
        'MarkerEdgeColor', [1 0 0], 'MarkerFaceColor', [1 0 0]);
    title('Velocity (peaks)')
    xlabel('contour');
    ylabel('time');
    xlim([0, prm]);
    ylim([0, N_frames]);
    
    subplot(1,4,4);
    hold on;
    axis ij;
    imagesc(C_tr);
    scatter(J_peak, I_peak, ...
        'Marker', 'o', 'SizeData', 5, ...
        'MarkerEdgeColor', [1 0 0], 'MarkerFaceColor', [1 0 0]);
    title('Biosensor Signal (peaks)');
    xlabel('contour');
    ylabel('time');
    xlim([0, prm]);
    ylim([0, N_frames]);
    saveas(fig, fullfile(save_fold, 'peaks_kymographs.png'));
    close(fig);
    
    fig = figure('Position', [50 50 1800 500]);
    subplot(1,2,1);
    hold on;
    grid on;
    box on;
    for i = 1:size(vw_peak_vals,1)
        plot(-w:1:w, vw_peak_vals(i,:), 'Color', [1 0 0 0.2], 'LineWidth', 1, 'LineStyle', '-', 'Marker', 'none');
    end
    plot(-w:1:w, mean(vw_peak_vals,1), 'Color', [0 0 0 1], 'LineWidth', 3, 'LineStyle', '-', 'Marker', 'none');
    plot(-w:1:w, zeros(1,2*w+1), 'Color', [0 1 0 1], 'LineWidth', 1, 'LineStyle', '--', 'Marker', 'none');
    title('velocity near peaks');
    
    subplot(1,2,2);
    hold on;
    grid on;
    box on;
    for i = 1:size(vw_peak_vals,1)
        plot(-w:1:w, Cw_peak_vals(i,:), 'Color', [1 0 0 0.1], 'LineWidth', 1, 'LineStyle', '-', 'Marker', 'none');
    end
    plot(-w:1:w, mean(Cw_peak_vals,1), 'Color', [0 0 0 1], 'LineWidth', 3, 'LineStyle', '-', 'Marker', 'none');
    title('Biosensor Signal near peaks');
    xlabel('frame');
    saveas(fig, fullfile(save_fold, 'peaks_velocity_concentration.png'));
    close(fig);
    
    fig = figure('Position', [50 50 1800 500]);
    subplot(1,4,1);
    hold on;
    axis ij;
    imagesc(sVEL);
    title('Velocity (smoothed)')
    xlabel('contour');
    ylabel('time');
    xlim([0, prm]);
    ylim([0, N_frames]);
    
    subplot(1,4,2);
    hold on;
    axis ij;
    imagesc(sVEL_dip_th);
    title('Velocity threshold')
    xlabel('contour');
    ylabel('time');
    xlim([0, prm]);
    ylim([0, N_frames]);
    
    subplot(1,4,3);
    hold on;
    axis ij;
    imagesc(v_tr);
    scatter(J_dip, I_dip, ...
        'Marker', 'o', 'SizeData', 5, ...
        'MarkerEdgeColor', [1 0 0], 'MarkerFaceColor', [1 0 0]);
    title('Velocity (dips)')
    xlabel('contour');
    ylabel('time');
    xlim([0, prm]);
    ylim([0, N_frames]);
    
    subplot(1,4,4);
    hold on;
    axis ij;
    imagesc(C_tr);
    scatter(J_dip, I_dip, ...
        'Marker', 'o', 'SizeData', 5, ...
        'MarkerEdgeColor', [1 0 0], 'MarkerFaceColor', [1 0 0]);
    title('Biosensor Signal (dips)');
    xlabel('contour');
    ylabel('time');
    xlim([0, prm]);
    ylim([0, N_frames]);
    saveas(fig, fullfile(save_fold, 'dips_kymographs.png'));
    close(fig);
    
    fig = figure('Position', [50 50 1800 500]);
    subplot(1,2,1);
    hold on;
    grid on;
    box on;
    for i = 1:size(vw_dip_vals,1)
        plot(-w:1:w, vw_dip_vals(i,:), 'Color', [1 0 0 0.1], 'LineWidth', 1, 'LineStyle', '-', 'Marker', 'none');
    end
    plot(-w:1:w, mean(vw_dip_vals,1), 'Color', [0 0 0 1], 'LineWidth', 3, 'LineStyle', '-', 'Marker', 'none');
    plot(-w:1:w, zeros(1,2*w+1), 'Color', [0 1 0 1], 'LineWidth', 1, 'LineStyle', '--', 'Marker', 'none');
    title('velocity near dips');
    
    subplot(1,2,2);
    hold on;
    grid on;
    box on;
    for i = 1:size(vw_dip_vals,1)
        plot(-w:1:w, Cw_dip_vals(i,:), 'Color', [1 0 0 0.1], 'LineWidth', 1, 'LineStyle', '-', 'Marker', 'none');
    end
    plot(-w:1:w, mean(Cw_dip_vals,1), 'Color', [0 0 0 1], 'LineWidth', 3, 'LineStyle', '-', 'Marker', 'none');
    title('Biosensor Signal near dips');
    xlabel('frame');
    saveas(fig, fullfile(save_fold, 'dips_velocity_concentration.png'));
    close(fig);
    
    fig = figure('Position', [50 50 800 800]);
    hold on;
    grid on;
    box on;
    
    colororder({'k','r'});
    yyaxis left;
    p1 = plot(-w:1:w, mean(vw_peak_vals,1), 'LineWidth', 3, 'LineStyle', '-', 'Marker', 'none', 'Color', [0 0 0 0.5]);
    plot(-w:1:w, vw_peak_CI_1, 'LineWidth', 1, 'LineStyle', '--', 'Marker', 'none', 'Color', [0 0 0 0.5]);
    plot(-w:1:w, vw_peak_CI_2, 'LineWidth', 1, 'LineStyle', '--', 'Marker', 'none', 'Color', [0 0 0 0.5]);
    x = [-w:1:w, fliplr(-w:1:w)];
    y = [vw_peak_CI_1, fliplr(vw_peak_CI_2)];
    patch(x, y, 'g', 'EdgeAlpha', 0.3, 'EdgeColor', [0 0 0], 'FaceAlpha', 0.3, 'FaceColor', [0 0 0]);
    xlim([-20, 20]);
    
    yyaxis right;
    p2 = plot(-w:1:w, mean(Cw_peak_vals,1), 'LineWidth', 3, 'LineStyle', '-', 'Marker', 'none', 'Color', [1 0 0 0.5]);
    plot(-w:1:w, Cw_peak_CI_1, 'LineWidth', 1, 'LineStyle', '--', 'Marker', 'none', 'Color', [1 0 0 0.5]);
    plot(-w:1:w, Cw_peak_CI_2, 'LineWidth', 1, 'LineStyle', '--', 'Marker', 'none', 'Color', [1 0 0 0.5]);
    x = [-w:1:w, fliplr(-w:1:w)];
    y = [Cw_peak_CI_1, fliplr(Cw_peak_CI_2)];
    patch(x, y, 'g', 'EdgeAlpha', 0.3, 'EdgeColor', [1 0 0], 'FaceAlpha', 0.3, 'FaceColor', [1 0 0]);
    xlim([-20, 20]);
    
    xlabel('Frame');
    legend([p1, p2], {'velocity', 'biosensor signal'});
    title('Mean velocity and biosensor signal near peaks');
    xlim([-w, w]);
    set(gca, 'FontSize', 15);
    saveas(fig, fullfile(save_fold, 'peaks_CI_velocity_concentration.png'));
    close(fig);
    
    fig = figure('Position', [50 50 800 800]);
    hold on;
    grid on;
    box on;
    
    colororder({'k','r'});
    yyaxis left;
    p1 = plot(-w:1:w, mean(vw_dip_vals,1), 'LineWidth', 3, 'LineStyle', '-', 'Marker', 'none', 'Color', [0 0 0 0.5]);
    plot(-w:1:w, vw_dip_CI_1, 'LineWidth', 1, 'LineStyle', '--', 'Marker', 'none', 'Color', [0 0 0 0.5]);
    plot(-w:1:w, vw_dip_CI_2, 'LineWidth', 1, 'LineStyle', '--', 'Marker', 'none', 'Color', [0 0 0 0.5]);
    x = [-w:1:w, fliplr(-w:1:w)];
    y = [vw_dip_CI_1, fliplr(vw_dip_CI_2)];
    patch(x, y, 'g', 'EdgeAlpha', 0.3, 'EdgeColor', [0 0 0], 'FaceAlpha', 0.3, 'FaceColor', [0 0 0]);
    xlim([-20, 20]);
    
    yyaxis right;
    p2 = plot(-w:1:w, mean(Cw_dip_vals,1), 'LineWidth', 3, 'LineStyle', '-', 'Marker', 'none', 'Color', [1 0 0 0.5]);
    plot(-w:1:w, Cw_dip_CI_1, 'LineWidth', 1, 'LineStyle', '--', 'Marker', 'none', 'Color', [1 0 0 0.5]);
    plot(-w:1:w, Cw_dip_CI_2, 'LineWidth', 1, 'LineStyle', '--', 'Marker', 'none', 'Color', [1 0 0 0.5]);
    x = [-w:1:w, fliplr(-w:1:w)];
    y = [Cw_dip_CI_1, fliplr(Cw_dip_CI_2)];
    patch(x, y, 'g', 'EdgeAlpha', 0.3, 'EdgeColor', [1 0 0], 'FaceAlpha', 0.3, 'FaceColor', [1 0 0]);
    xlim([-20, 20]);
    
    xlabel('Frame');
    legend([p1, p2], {'velocity', 'biosensor signal'});
    title('Mean velocity and biosensor signal near dips');
    xlim([-w, w]);
    set(gca, 'FontSize', 15);
    saveas(fig, fullfile(save_fold, 'dips_CI_velocity_concentration.png'));
    close(fig);
    
    
    fig = figure('Position', [50 50 800 800]);
    hold on;
    grid on;
    box on;
    
    colororder({'k','r'});
    yyaxis left;
    p1 = plot(-w:1:w, mean(xw_peak_vals,1), 'LineWidth', 3, 'LineStyle', '-', 'Marker', 'none', 'Color', [0 0 0 0.5]);
    plot(-w:1:w, xw_peak_CI_1, 'LineWidth', 1, 'LineStyle', '--', 'Marker', 'none', 'Color', [0 0 0 0.5]);
    plot(-w:1:w, xw_peak_CI_2, 'LineWidth', 1, 'LineStyle', '--', 'Marker', 'none', 'Color', [0 0 0 0.5]);
    x = [-w:1:w, fliplr(-w:1:w)];
    y = [xw_peak_CI_1, fliplr(xw_peak_CI_2)];
    patch(x, y, 'g', 'EdgeAlpha', 0.3, 'EdgeColor', [0 0 0], 'FaceAlpha', 0.3, 'FaceColor', [0 0 0]);
    xlim([-20, 20]);
    
    yyaxis right;
    p2 = plot(-w:1:w, mean(Cw_peak_vals,1), 'LineWidth', 3, 'LineStyle', '-', 'Marker', 'none', 'Color', [1 0 0 0.5]);
    plot(-w:1:w, Cw_peak_CI_1, 'LineWidth', 1, 'LineStyle', '--', 'Marker', 'none', 'Color', [1 0 0 0.5]);
    plot(-w:1:w, Cw_peak_CI_2, 'LineWidth', 1, 'LineStyle', '--', 'Marker', 'none', 'Color', [1 0 0 0.5]);
    x = [-w:1:w, fliplr(-w:1:w)];
    y = [Cw_peak_CI_1, fliplr(Cw_peak_CI_2)];
    patch(x, y, 'g', 'EdgeAlpha', 0.3, 'EdgeColor', [1 0 0], 'FaceAlpha', 0.3, 'FaceColor', [1 0 0]);
    xlim([-20, 20]);
    
    xlabel('Frame');
    legend([p1, p2], {'coordinate', 'biosensor signal'});
    title('Mean coordinate and biosensor signal near peaks');
    xlim([-w, w]);
    set(gca, 'FontSize', 15);
    saveas(fig, fullfile(save_fold, 'peaks_CI_coordinate_concentration.png'));
    close(fig);
    
    
    fig = figure('Position', [50 50 800 800]);
    hold on;
    grid on;
    box on;
    
    colororder({'k','r'});
    yyaxis left;
    p1 = plot(-w:1:w, mean(xw_dip_vals,1), 'LineWidth', 3, 'LineStyle', '-', 'Marker', 'none', 'Color', [0 0 0 0.5]);
    plot(-w:1:w, xw_dip_CI_1, 'LineWidth', 1, 'LineStyle', '--', 'Marker', 'none', 'Color', [0 0 0 0.5]);
    plot(-w:1:w, xw_dip_CI_2, 'LineWidth', 1, 'LineStyle', '--', 'Marker', 'none', 'Color', [0 0 0 0.5]);
    x = [-w:1:w, fliplr(-w:1:w)];
    y = [xw_dip_CI_1, fliplr(xw_dip_CI_2)];
    patch(x, y, 'g', 'EdgeAlpha', 0.3, 'EdgeColor', [0 0 0], 'FaceAlpha', 0.3, 'FaceColor', [0 0 0]);
    xlim([-20, 20]);
    
    yyaxis right;
    p2 = plot(-w:1:w, mean(Cw_dip_vals,1), 'LineWidth', 3, 'LineStyle', '-', 'Marker', 'none', 'Color', [1 0 0 0.5]);
    plot(-w:1:w, Cw_dip_CI_1, 'LineWidth', 1, 'LineStyle', '--', 'Marker', 'none', 'Color', [1 0 0 0.5]);
    plot(-w:1:w, Cw_dip_CI_2, 'LineWidth', 1, 'LineStyle', '--', 'Marker', 'none', 'Color', [1 0 0 0.5]);
    x = [-w:1:w, fliplr(-w:1:w)];
    y = [Cw_dip_CI_1, fliplr(Cw_dip_CI_2)];
    patch(x, y, 'g', 'EdgeAlpha', 0.3, 'EdgeColor', [1 0 0], 'FaceAlpha', 0.3, 'FaceColor', [1 0 0]);
    xlim([-20, 20]);
    
    xlabel('Frame');
    legend([p1, p2], {'coordinate', 'biosensor signal'});
    title('Mean velocity and biosensor signal near dips');
    xlim([-w, w]);
    set(gca, 'FontSize', 15);
    saveas(fig, fullfile(save_fold, 'dips_CI_coordinate_concentration.png'));
    close(fig);
    
    
    fig = figure('Position', [50 50 800 800]);
    hold on;
    grid on;
    box on;
    
    colororder({'k','r'});
    yyaxis left;
    p1 = plot(-w+diff_w:1:(w-diff_w), mean(vw_peak_vals(:,1+diff_w:end-diff_w),1), 'LineWidth', 3, 'LineStyle', '-', 'Marker', 'none', 'Color', [0 0 0 0.5]);
    plot(-w+diff_w:1:(w-diff_w), vw_peak_CI_1(1+diff_w:end-diff_w), 'LineWidth', 1, 'LineStyle', '--', 'Marker', 'none', 'Color', [0 0 0 0.5]);
    plot(-w+diff_w:1:(w-diff_w), vw_peak_CI_2(1+diff_w:end-diff_w), 'LineWidth', 1, 'LineStyle', '--', 'Marker', 'none', 'Color', [0 0 0 0.5]);
    x = [-w+diff_w:1:(w-diff_w), fliplr(-w+diff_w:1:(w-diff_w))];
    y = [vw_peak_CI_1(1+diff_w:end-diff_w), fliplr(vw_peak_CI_2(1+diff_w:end-diff_w))];
    patch(x, y, 'g', 'EdgeAlpha', 0.3, 'EdgeColor', [0 0 0], 'FaceAlpha', 0.3, 'FaceColor', [0 0 0]);
    xlim([-w, (w-diff_w)]);
    
    yyaxis right;
    p2 = plot(-w+diff_w:1:(w-diff_w), mean(dCw_peak_vals,1), 'LineWidth', 3, 'LineStyle', '-', 'Marker', 'none', 'Color', [1 0 0 0.5]);
    plot(-w+diff_w:1:(w-diff_w), dCw_peak_CI_1, 'LineWidth', 1, 'LineStyle', '--', 'Marker', 'none', 'Color', [1 0 0 0.5]);
    plot(-w+diff_w:1:(w-diff_w), dCw_peak_CI_2, 'LineWidth', 1, 'LineStyle', '--', 'Marker', 'none', 'Color', [1 0 0 0.5]);
    x = [-w+diff_w:1:(w-diff_w), fliplr(-w+diff_w:1:(w-diff_w))];
    y = [dCw_peak_CI_1, fliplr(dCw_peak_CI_2)];
    patch(x, y, 'g', 'EdgeAlpha', 0.3, 'EdgeColor', [1 0 0], 'FaceAlpha', 0.3, 'FaceColor', [1 0 0]);
    xlim([-w+diff_w, (w-diff_w)]);
    
    xlabel('Frame');
    legend([p1, p2], {'velocity', 'diff(biosensor signal)'});
    title('Mean velocity and diff(biosensor signal) near peaks');
    set(gca, 'FontSize', 15);

    saveas(fig, fullfile(save_fold, 'peaks_CI_velocity_diff_concentration.png'));
    close(fig);
    
    fig = figure('Position', [50 50 800 800]);
    hold on;
    grid on;
    box on;
    
    colororder({'k','r'});
    yyaxis left;
    p1 = plot(-w+diff_w:1:(w-diff_w), mean(vw_dip_vals(:,1+diff_w:end-diff_w),1), 'LineWidth', 3, 'LineStyle', '-', 'Marker', 'none', 'Color', [0 0 0 0.5]);
    plot(-w+diff_w:1:(w-diff_w), vw_dip_CI_1(1+diff_w:end-diff_w), 'LineWidth', 1, 'LineStyle', '--', 'Marker', 'none', 'Color', [0 0 0 0.5]);
    plot(-w+diff_w:1:(w-diff_w), vw_dip_CI_2(1+diff_w:end-diff_w), 'LineWidth', 1, 'LineStyle', '--', 'Marker', 'none', 'Color', [0 0 0 0.5]);
    x = [-w+diff_w:1:(w-diff_w), fliplr(-w+diff_w:1:(w-diff_w))];
    y = [vw_dip_CI_1(1+diff_w:end-diff_w), fliplr(vw_dip_CI_2(1+diff_w:end-diff_w))];
    patch(x, y, 'g', 'EdgeAlpha', 0.3, 'EdgeColor', [0 0 0], 'FaceAlpha', 0.3, 'FaceColor', [0 0 0]);
    xlim([-w, (w-diff_w)]);
    
    yyaxis right;
    p2 = plot(-w+diff_w:1:(w-diff_w), mean(dCw_dip_vals,1), 'LineWidth', 3, 'LineStyle', '-', 'Marker', 'none', 'Color', [1 0 0 0.5]);
    plot(-w+diff_w:1:(w-diff_w), dCw_dip_CI_1, 'LineWidth', 1, 'LineStyle', '--', 'Marker', 'none', 'Color', [1 0 0 0.5]);
    plot(-w+diff_w:1:(w-diff_w), dCw_dip_CI_2, 'LineWidth', 1, 'LineStyle', '--', 'Marker', 'none', 'Color', [1 0 0 0.5]);
    x = [-w+diff_w:1:(w-diff_w), fliplr(-w+diff_w:1:(w-diff_w))];
    y = [dCw_dip_CI_1, fliplr(dCw_dip_CI_2)];
    patch(x, y, 'g', 'EdgeAlpha', 0.3, 'EdgeColor', [1 0 0], 'FaceAlpha', 0.3, 'FaceColor', [1 0 0]);
    xlim([-w+diff_w, (w-diff_w)]);
    
    xlabel('Frame');
    legend([p1, p2], {'velocity', 'diff(biosensor signal)'});
    title('Mean velocity and diff(biosensor signal) near dips');
    set(gca, 'FontSize', 15);

    saveas(fig, fullfile(save_fold, 'dips_CI_velocity_diff_concentration.png'));
    close(fig);

    fig = figure('Position', [50 50 800 800]);
    hold on;
    grid on;
    box on;
    for p = 1:prm
        plot(-w:1:w, corr_vel_C(:,p), 'Color', [1 0 0 0.05], 'LineWidth', 1, 'LineStyle', '-', 'Marker', 'none');
    end
    plot(-w:1:w, mean(corr_vel_C,2), 'Color', [0 0 0 1], 'LineWidth', 3, 'LineStyle', '-', 'Marker', 'none');
    plot(-w:1:w, zeros(1,2*w+1), 'Color', [0 1 0 1], 'LineWidth', 1, 'LineStyle', '--', 'Marker', 'none');
    plot(zeros(1,2*w+1), linspace(-0.5, 0.5, 2*w+1), 'Color', [0 1 0 1], 'LineWidth', 1, 'LineStyle', '--', 'Marker', 'none');
    title('Correlation velocity vs. Biosensor Signal');
    xlabel('Time lag (frames)');
    saveas(fig, fullfile(save_fold, 'crosscorrelations_velocity_vs_concentration.png'));
    close(fig);
    
%     fig = figure('Position', [50 50 500 500]);
%     hold on;
%     grid on;
%     box on;
%     colororder({'r','g'});
%     xlim([min(v_nodes), max(v_nodes)]);
%     xticks(v_nodes);
%     set(gca, 'FontSize', 15);
%     xlabel('Velocity');
%     ylabel('Biosensor signal', 'Color', [0 0 0]);
%     
%     plot(v_nodes, C_tr_all_mean_vals, 'LineWidth', 2, 'LineStyle', '-', 'Marker', 'none', 'Color', [1 0 0 0.5]);
%     plot(v_nodes, C_tr_all_CI_1_vals, 'LineWidth', 1, 'LineStyle', '--', 'Marker', 'none', 'Color', [1 0 0 0.5]);
%     plot(v_nodes, C_tr_all_CI_2_vals, 'LineWidth', 1, 'LineStyle', '--', 'Marker', 'none', 'Color', [1 0 0 0.5]);
%     x = [v_nodes, fliplr(v_nodes)];
%     y = [C_tr_all_CI_1_vals, fliplr(C_tr_all_CI_2_vals)];
%     patch(x, y, 'r', 'EdgeAlpha', 0.3, 'EdgeColor', [1 0 0], 'FaceAlpha', 0.3, 'FaceColor', [1 0 0]);
%     title('Velocity vs. Biosensor Signal (all)');
%     drawnow;
%     saveas(fig, fullfile(save_fold, 'velocity_vs_concentration_averaged_all_trajectories.png'));
%     close(fig);
%     
%     fig = figure('Position', [50 50 500 500]);
%     hold on;
%     grid on;
%     box on;
%     colororder({'r','g'});
%     xlim([min(v_nodes), max(v_nodes)]);
%     xticks(v_nodes);
%     set(gca, 'FontSize', 15);
%     xlabel('Velocity');
%     ylabel('Biosensor signal', 'Color', [0 0 0]);
%     
%     plot(v_nodes, C_tr_peak_mean_vals, 'LineWidth', 2, 'LineStyle', '-', 'Marker', 'none', 'Color', [1 0 0 0.5]);
%     plot(v_nodes, C_tr_peak_CI_1_vals, 'LineWidth', 1, 'LineStyle', '--', 'Marker', 'none', 'Color', [1 0 0 0.5]);
%     plot(v_nodes, C_tr_peak_CI_2_vals, 'LineWidth', 1, 'LineStyle', '--', 'Marker', 'none', 'Color', [1 0 0 0.5]);
%     x = [v_nodes, fliplr(v_nodes)];
%     y = [C_tr_peak_CI_1_vals, fliplr(C_tr_peak_CI_2_vals)];
%     patch(x, y, 'r', 'EdgeAlpha', 0.3, 'EdgeColor', [1 0 0], 'FaceAlpha', 0.3, 'FaceColor', [1 0 0]);
%     title('Velocity vs. Biosensor Signal (peaks)');
%     drawnow;
%     saveas(fig, fullfile(save_fold, 'velocity_vs_concentration_averaged_peaks.png'));
%     close(fig);
%     
%     fig = figure('Position', [50 50 500 500]);
%     hold on;
%     grid on;
%     box on;
%     colororder({'r','g'});
%     xlim([min(v_nodes), max(v_nodes)]);
%     xticks(v_nodes);
%     set(gca, 'FontSize', 15);
%     xlabel('Velocity');
%     ylabel('Biosensor signal', 'Color', [0 0 0]);
%     
%     plot(v_nodes, C_tr_dip_mean_vals, 'LineWidth', 2, 'LineStyle', '-', 'Marker', 'none', 'Color', [1 0 0 0.5]);
%     plot(v_nodes, C_tr_dip_CI_1_vals, 'LineWidth', 1, 'LineStyle', '--', 'Marker', 'none', 'Color', [1 0 0 0.5]);
%     plot(v_nodes, C_tr_dip_CI_2_vals, 'LineWidth', 1, 'LineStyle', '--', 'Marker', 'none', 'Color', [1 0 0 0.5]);
%     x = [v_nodes, fliplr(v_nodes)];
%     y = [C_tr_dip_CI_1_vals, fliplr(C_tr_dip_CI_2_vals)];
%     patch(x, y, 'r', 'EdgeAlpha', 0.3, 'EdgeColor', [1 0 0], 'FaceAlpha', 0.3, 'FaceColor', [1 0 0]);
%     title('Velocity vs. Biosensor Signal (dips)');
%     drawnow;
%     saveas(fig, fullfile(save_fold, 'velocity_vs_concentration_averaged_dips.png'));
%     close(fig);
%     
%     
%     [v_tr_all_vals_smpl, C_tr_all_vals_smpl, C_dens_vals] = compute_density(v_tr_all_vals, C_tr_all_vals, v_dw, C_dw, rand_frac);
%     c_data = color_data(C_dens_vals, cm);
%     
%     fig = figure('Position', [50 50 500 500]);
%     hold on;
%     grid on;
%     box on;
%     scatter(v_tr_all_vals_smpl, C_tr_all_vals_smpl, 'Marker', '.', 'CData', c_data);
%     xlabel('Velocity');
%     ylabel('Biosensor signal');
%     set(gca, 'FontSize', 15);
%     xlim([min(v_nodes), max(v_nodes)]);
%     title('All trajectories');
%     drawnow;
%     saveas(fig, fullfile(save_fold, 'velocity_vs_concentration_cloud_all_trajectories.png'));
%     close(fig);
%     
%     [v_tr_peak_vals_smpl, C_tr_peak_vals_smpl, C_dens_vals] = compute_density(v_tr_peak_vals, C_tr_peak_vals, v_dw, C_dw, rand_frac);
%     c_data = color_data(C_dens_vals, cm);
%     
%     fig = figure('Position', [50 50 500 500]);
%     hold on;
%     grid on;
%     box on;
%     scatter(v_tr_peak_vals_smpl, C_tr_peak_vals_smpl, 'Marker', '.', 'CData', c_data);
%     xlabel('Velocity');
%     ylabel('Biosensor signal');
%     set(gca, 'FontSize', 15);
%     xlim([min(v_nodes), max(v_nodes)]);
%     title('Peaks');
%     drawnow;
%     saveas(fig, fullfile(save_fold, 'velocity_vs_concentration_cloud_peaks.png'));
%     close(fig);
%     
%     [v_tr_dip_vals_smpl, C_tr_dip_vals_smpl, C_dens_vals] = compute_density(v_tr_dip_vals, C_tr_dip_vals, v_dw, C_dw, rand_frac);
%     c_data = color_data(C_dens_vals, cm);
%     
%     fig = figure('Position', [50 50 500 500]);
%     hold on;
%     grid on;
%     box on;
%     scatter(v_tr_dip_vals_smpl, C_tr_dip_vals_smpl, 'Marker', '.', 'CData', c_data);
%     xlabel('Velocity');
%     ylabel('Biosensor signal');
%     set(gca, 'FontSize', 15);
%     xlim([min(v_nodes), max(v_nodes)]);
%     title('Dips');
%     drawnow;
%     saveas(fig, fullfile(save_fold, 'velocity_vs_concentration_cloud_dips.png'));
%     close(fig);
%     
%     
%     ids = v_tr_all_vals>0;
%     v_tr_all_pos_vals = v_tr_all_vals(ids);
%     C_tr_all_pos_vals = C_tr_all_vals(ids);
%     
%     ids = v_tr_all_vals<0;
%     v_tr_all_neg_vals = v_tr_all_vals(ids);
%     C_tr_all_neg_vals = C_tr_all_vals(ids);
%     
%     lm_all_pos = fitlm(v_tr_all_pos_vals, C_tr_all_pos_vals);
%     lm_all_neg = fitlm(v_tr_all_neg_vals, C_tr_all_neg_vals);
%     
%     fig = figure('Position', [50 50 500 500]);
%     hold on;
%     grid on;
%     box on;
%     pl = plot(lm_all_pos);
%     set(pl(2), 'LineWidth', 3);
%     set(pl(3), 'LineWidth', 2);
%     set(pl(4), 'LineWidth', 2);
%     pl = plot(lm_all_neg);
%     set(pl(2), 'LineWidth', 3);
%     set(pl(3), 'LineWidth', 2);
%     set(pl(4), 'LineWidth', 2);
%     title('All trajectories');
%     xlabel('Velocity');
%     ylabel('Biosensor signal');
%     set(gca, 'FontSize', 15);
%     drawnow;
%     saveas(fig, fullfile(save_fold, 'velocity_vs_concentration_cloud_linear_regression_all_trajectories.png'));
%     close(fig);
%     
%     ids = v_tr_peak_vals>0;
%     v_tr_peak_pos_vals = v_tr_peak_vals(ids);
%     C_tr_peak_pos_vals = C_tr_peak_vals(ids);
%     
%     ids = v_tr_peak_vals<0;
%     v_tr_peak_neg_vals = v_tr_peak_vals(ids);
%     C_tr_peak_neg_vals = C_tr_peak_vals(ids);
%     
%     lm_peak_pos = fitlm(v_tr_peak_pos_vals, C_tr_peak_pos_vals);
%     lm_peak_neg = fitlm(v_tr_peak_neg_vals, C_tr_peak_neg_vals);
%     
%     fig = figure('Position', [50 50 500 500]);
%     hold on;
%     grid on;
%     box on;
%     pl = plot(lm_peak_pos);
%     set(pl(2), 'LineWidth', 3);
%     set(pl(3), 'LineWidth', 2);
%     set(pl(4), 'LineWidth', 2);
%     pl = plot(lm_peak_neg);
%     set(pl(2), 'LineWidth', 3);
%     set(pl(3), 'LineWidth', 2);
%     set(pl(4), 'LineWidth', 2);
%     title('Peaks');
%     xlabel('Velocity');
%     ylabel('Biosensor signal');
%     set(gca, 'FontSize', 15);
%     drawnow;
%     saveas(fig, fullfile(save_fold, 'velocity_vs_concentration_cloud_linear_regression_peaks.png'));
%     close(fig);
%     
%     ids = v_tr_dip_vals>0;
%     v_tr_dip_pos_vals = v_tr_dip_vals(ids);
%     C_tr_dip_pos_vals = C_tr_dip_vals(ids);
%     
%     ids = v_tr_dip_vals<0;
%     v_tr_dip_neg_vals = v_tr_dip_vals(ids);
%     C_tr_dip_neg_vals = C_tr_dip_vals(ids);
%     
%     lm_dip_pos = fitlm(v_tr_dip_pos_vals, C_tr_dip_pos_vals);
%     lm_dip_neg = fitlm(v_tr_dip_neg_vals, C_tr_dip_neg_vals);
%     
%     fig = figure('Position', [50 50 500 500]);
%     hold on;
%     grid on;
%     box on;
%     pl = plot(lm_dip_pos);
%     set(pl(2), 'LineWidth', 3);
%     set(pl(3), 'LineWidth', 2);
%     set(pl(4), 'LineWidth', 2);
%     pl = plot(lm_dip_neg);
%     set(pl(2), 'LineWidth', 3);
%     set(pl(3), 'LineWidth', 2);
%     set(pl(4), 'LineWidth', 2);
%     title('Dips');
%     xlabel('Velocity');
%     ylabel('Biosensor signal');
%     set(gca, 'FontSize', 15);
%     drawnow;
%     saveas(fig, fullfile(save_fold, 'velocity_vs_concentration_cloud_linear_regression_dips.png'));
%     close(fig);
    
    save(fullfile(save_fold, 'inspect_trajectory_data.mat'), ...
        'vw_peak_vals', 'Cw_peak_vals', 'xw_peak_vals', 'dCw_peak_vals', ...
        'vw_dip_vals', 'Cw_dip_vals', 'xw_dip_vals', 'dCw_dip_vals', ...
        'corr_vel_C', 'v_vals', 'C_vals', ...
        'v_vals2', 'C_vals2', 'C_stds2', 'diff_w');
end

