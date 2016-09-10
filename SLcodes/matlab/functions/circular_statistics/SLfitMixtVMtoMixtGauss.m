

%SLfitMixtVMtoMixtGauss.m
%
% author: Steeve Laquitaine
%   date: 150123
%purpose: least square fit a miture of von mises to a mixture of Gaussian
%         distribution
%
%usage:
%
%  [SSE,strength,exitflag,outputFit] = SLfitMixtVMtoMixtGauss(12,[145 305],107,5:10:355,100)


function [SSE,strength,exitflag,outputFit] = SLfitMixtVMtoMixtGauss(strength,modes,trialnum,x,str0)

%Mixture of Gaussian
y = SLmakeDiscreteMixtureGauss([strength strength],5:10:355,[145 305],107);
yMoG = y.parameter.dir.count;

%least square fit 
mode1 = modes(1);
mode2 = modes(2);

%fit
options = optimset('TolFun',0);

%least square fit 
[strength,SSE,exitflag,outputFit] = fminsearch(@(strength) ...
    makeSSE(x,yMoG,strength,mode1,mode2,trialnum),...
    str0,options);

%Best predictions
MoVM = makeDiscreteMixtureVM(x,strength,mode1,mode2,trialnum);

%draw best predictions
clf
set(gcf,'color','w')
hold all;
plot(x,yMoG,'-ro')
plot(x,MoVM,'--ko')
legend('MoG','MoVM')
legend('boxoff')


%Sum of squared error
function [SSE,strength] = makeSSE(x,yMoG,strength,mode1,mode2,trialnum)

%predited
MoVM = makeDiscreteMixtureVM(x,strength,mode1,mode2,trialnum);

%draw
set(gcf,'color','w')
hold all;
plot(x,MoVM,'-ko')
plot(x,yMoG,'-ro')
legend('MoVM','MoG')
legend('boxoff')
drawnow

%SSE
SSE = sum((yMoG - MoVM).^2);

%update
fprintf('%.2f  %.2f \n',SSE,strength)


%make Discrete Mixture von Mises
function count = makeDiscreteMixtureVM(x,strength,mode1,mode2,trialnum)

%initialize
series = [];
count  = [];
sampsiz = numel(x);

%Bimodal prior
if strength~=0
    
    %Mixture of two von Mises
    MoVM = SLMixtureVM(x,[mode1 mode2],strength); 
    
    %count of each direction
    f2 = MoVM * trialnum ;
    
    %Repeat and store samples
    for i = 1 : numel(x)
        f2repet = [];
        f2repet = repmat(x(i),round(f2(i)),1);
        series = [series; f2repet];
    end
    series = series';
end

%count
for i = 1 : sampsiz    
    count(i) = numel(find(series==x(i)));
end
