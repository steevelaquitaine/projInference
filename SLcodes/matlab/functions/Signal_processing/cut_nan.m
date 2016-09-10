function y = cut_nan(x)

y=x(1:max(max(size(x,1)-sum(isnan(x),1))),:)