function p=test_smooth(x)

p=nan(size(x,1),size(x,2));
sigma=5%  default: 3;

col_size=size(x,2);
row_size=size(x,1);

for i2=1:col_size; % col
      for i1=1:row_size %line
          range_filter=row_size;
          filter_tmp=(-range_filter:range_filter)';
          filter=exp(-filter_tmp.^2/(2*sigma^2));
          p(i1,i2)=nansum(x(:,i2).*filter(range_filter-i1+2:2*range_filter-i1+1)) ./nansum(filter(range_filter-i1+2:2*range_filter-i1+1)); 
     end
end
p(find(isnan(x)==1))=NaN;