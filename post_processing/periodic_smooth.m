function [result] = periodic_smooth(data,window_len)
    data1 = [data(:); data(:); data(:)];
    data2 = smooth(data1,window_len);
    result = zeros(size(data));
    l = length(data);
    result(:) = data2(l+1:l+l);
end