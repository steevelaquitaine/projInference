function p=m_nanngmaw_col(x)
% steeve last version 12/ 10/2008 

% use a matrix
% smooth the matrix columns with a gaussian moving window normalised over the vector length

clc;


%x=cut_nan(x);
p=nan(size(x,1),size(x,2));

sigma=1.2%  default: 3;
for i2=1:length(x(1,:)) % col
    x1(i2)=length(x)%(find(isnan(x(:,i2))==0)); 
    for i1=1:x1(i2); % line 
    range_filter= x1(i2);
    filter_tmp=(-range_filter:range_filter)';
    filter=exp(-filter_tmp.^2/(2*sigma^2));
    p(i1,i2)=nansum(x(1:x1(i2),i2).*filter(range_filter-i1+2:2*range_filter-i1+1))./nansum(filter(range_filter-i1+2:2*range_filter-i1+1));
    end
end

p(find(isnan(x)==1))=NaN;
clear('sigma','range_filter','filter_tmp','filter','i1','i2');


