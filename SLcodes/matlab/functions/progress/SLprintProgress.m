
 %author: Steeve Laquitaine
   %date: 140425
%purpose: print progress in percent
%usage: 
%for i=1:1000
%     SLprintProgress(i,1000)
%end

function SLprintProgress(i,numi)

%print progress
fprintf(1,'\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%s',['Progress ',num2str(round(i/numi*100)),'%']);
