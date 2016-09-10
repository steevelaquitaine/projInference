

% smooth data with a moving window you can chose
function p=plotstbin2(x)


sigma=1.3; 

length_filter=3
       
for i2=1:length(x(1,:)) % col
    x1(i2)=length(find(isnan(x(1:length_filter,i2))==0)); 
    for i1=1:length(x(:,1)); % line  
    length_filter=x1(i2); 
    
    filter_tmp=[-(length_filter-1)/2:(length_filter-1)/2];  
    filter=exp(-filter_tmp.^2/(2*sigma^2));                
    filter=filter/sum(filter); 
    
    x_for_filter(i1,:)=x(i1,end)*ones(1,length(x(i1,:))+length_filter-1);     
    x_for_filter(i1,1:(length_filter-1)/2)=x(i1,1);                
    x_for_filter(i1,(length_filter-1)/2+1:(length_filter-1)/2+length(x(1,:)))=x(i1,:); 
 
    p(i1,i2)=sum(filter.*x_for_filter(i1,[i2:i2+length_filter-1]));
    
    end
end



