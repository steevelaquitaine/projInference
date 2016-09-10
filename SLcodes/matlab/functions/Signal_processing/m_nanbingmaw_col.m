

% smooth data with a moving window you can chose
function p=nanbingmaw_col(z)

%nanclean
maxim=(sum(isnan(z)==0))';
[a,y]=max(maxim);
j=z(1:a,:);

%smooth
sigma=3; 
length_filter=10;
filter_tmp=[-(length_filter-1)/2:(length_filter-1)/2];  %-1:1
filter=(exp(-filter_tmp.^2/(2*sigma^2)))/sum(filter);  

for i1=1:length(j(:,1)); % line
    for i2=1:length(j(1,:)) % col 
    x1(i2)=length(find(isnan(j(1:length_filter,i2))==0));              
    j_for_filter(:,i2)=j(end,i2)*ones(length(j(:,i2))+length_filter-1,1); 
    j_for_filter(1:(length_filter-1)/2,i2)=j(1,i2); 
    j_for_filter((length_filter-1)/2+1:(length_filter-1)/2+length(j(:,1)),i2)=j(:,i2); 
    filter=filter';
    p(i1,i2)=sum(filter.*j_for_filter([i1:i1+length_filter-1],i2));  
    end
end

%initial
p(1:length_filter/2,:)=p(length_filter/2+1:length_filter,:);



