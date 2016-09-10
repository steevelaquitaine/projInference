
% plot Bined data function of many blocks of 'bin'trials
% data(:,1) = BinedData


% smooth data with a moving window you can chose
function m=bindata(z)


bin=50;
a=nan(size(z,1)/bin+1,1);
a(1,1)=1;
for i= 2:size(z,1)/bin+1
    a(i,1)=(i-1)*bin+1; 
end
m=nan(size(z,1)/bin,size(z,2))
for i=1:size(z,1)/bin
    m(i,:)=(mean(z(a(i):a(i+1),:)));
end

