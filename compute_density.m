function [v_vals, C_vals, dens_vals] = compute_density(v_vals, C_vals, v_dw, C_dw, rand_frac)
%COMPUTE_DENSITY Summary of this function goes here
%   Detailed explanation goes here

rand_ind = randsample(1:length(v_vals), fix(rand_frac*length(v_vals)));
v_vals = v_vals(rand_ind);
C_vals = C_vals(rand_ind);

dens_vals = zeros(length(v_vals), 1);
for i = 1:length(v_vals)
    v_i = v_vals(i);
    C_i = C_vals(i);
    ids = (v_vals > v_i - v_dw) & (v_vals < v_i + v_dw) & (C_vals > C_i - C_dw) & (C_vals < C_i + C_dw);
    dens_vals(i) = sum(ids);
end

end

