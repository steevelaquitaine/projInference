
 %author: steeve laquitaine
   %date: 140507
%purpose: draw softmax classification
  %usage:
    %x1=50;
    %x2=1:100;
    %sens=10;
    %PclassIsX2=SLdrawSoftmax(x1,x2,sens,Optiondisp);

function PclassIsX2=SLdrawSoftmax(x1,x2,sens,Optiondisp)

%iterate over values of x2 and use to calculate probability to classify as 
%x2. The greater x2's value the higher the probability to classify as x2.
PclassIsX2=nan(length(x2),1);
for j=1:length(x2)
    
    %Probability scaling factor
    scale(j)=exp(x2(j)/sens)+exp(x1/sens);
    
    %probability to be in class x1
    PclassIsX2(j)=exp(x2(j)/sens)./scale(j);
end

%draw
if strcmp(Optiondisp,'disp=on')
    figure('color','w')
    plot(x2-x1,PclassIsX2,'linesmoothing','on','linewidth',4,...
        'color',[.5 .5 .5])
    ylabel('P(class=X2)','fontsize',14)
    ylim([0 1])
    xlabel('x (decision variable)','fontsize',14)
    box off
    set(gca,'fontsize',14)
end
end