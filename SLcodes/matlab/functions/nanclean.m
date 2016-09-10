% smooth data with a moving window you can chose
function d=nanclean(A)

maxim=(sum(isnan(A)==0))';
[x,y]=max(maxim);
d=A(1:x,:);



