function c_data = color_data(dens_vals, cm)
%COLOR_DATA Summary of this function goes here
%   Detailed explanation goes here

c_N = size(cm, 1);
c_data = zeros(length(dens_vals), 3);
for i = 1:length(dens_vals)
    id = fix((dens_vals(i)-min(dens_vals))/(max(dens_vals)-min(dens_vals))*(c_N-1)) + 1;
    c_data(i,:) = cm(id,:);
end

end

