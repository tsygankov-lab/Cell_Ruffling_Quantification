function C_tr = read_concentration_kymograph(file_location, N_frames, prm, t_step, R_disk, x_tr, y_tr)
%READ_CONCENTRATION_KYMOGRAPH Summary of this function goes here
%   Detailed explanation goes here

C_tr = zeros(N_frames-t_step, prm);

se = strel('disk',R_disk,0);
nhood = se.Neighborhood;

for fr = 1:N_frames-t_step
    for t_step_i = 0:t_step
        x = x_tr(fr+t_step_i,:);
        y = y_tr(fr+t_step_i,:);
        C = double(imread(file_location, fr+t_step_i));
        xo = round(x);
        yo = round(y);
        for p = 1:prm
            if (yo(p)-R_disk) > 0 && (xo(p)-R_disk) > 0
                nhd = C((yo(p)-R_disk):(yo(p)+R_disk),(xo(p)-R_disk):(xo(p)+R_disk));
                nhd = nhd.*nhood;
                val = mean(nhd(nhd>0));
                if isnan(val)
                    val = 0;
                end
                C_tr(fr,p) = C_tr(fr,p) + val;
            end
        end
    end
end
C_tr(fr,p) = C_tr(fr,p)/(t_step+1);

end

