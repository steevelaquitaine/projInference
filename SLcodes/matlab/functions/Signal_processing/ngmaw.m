function [p]=ngmaw(x,sigma)
% steeve last version 12/ 10/2008 

% use a vector
% smooth the vector with a gaussian moving window normalised over the vector length


range_filter=length(x);
filter_tmp=(-range_filter:range_filter)';
filter=exp(-filter_tmp.^2/(2*sigma^2));

p=-ones(1,length(x));
for i1=1:length(x);
    p(i1)=nansum(x.*filter(range_filter-i1+2:2*range_filter-i1+1))./nansum(filter(range_filter-i1+2:2*range_filter-i1+1));
end

p=p';

clear('sigma','range_filter','filter_tmp','filter','i1');
