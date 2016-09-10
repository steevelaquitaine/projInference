


% rewrite 'repeat' times 'i*interval'

repeat=3
interval=10
for i=1:length(data)/repeat
    a=i*repeat; % start of the data vector (set) to average
    L((i-1)*repeat+1:i*repeat,1)=i*interval
end