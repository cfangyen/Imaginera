function vshift = linshift(v,n)
% % Right shift a 1-D array (assume to be a column vector)
a = zeros(size(v));
b = cat(1, a(1:n), v);
vshift = b(1:numel(v));

end

