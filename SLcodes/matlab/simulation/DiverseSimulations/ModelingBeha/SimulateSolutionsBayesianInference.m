function [ output_args ] = SimulateSolutionsBayesianInference( input_args )
%SIMULATEBAYESIANINFERENCE Summary of this function goes here
%   Detailed explanation goes here

%%
% for given ul and uP
ul = 5;
uP = 225;

% Loop over all pairs {sp, sl}
for sp = 1 : 50
    for sl = 1 : 50
        % Compute uPo and sPo
        uPo(sl,sp) = [(1./(1+(sl/sp).^2)).*ul + (1./(1+(sp/sl).^2)).*uP]';
        sPo(sl,sp) = sqrt(1/((1/sl)^2 + (1/sp)^2 ))';
    end
end

% Draw
screen = get(0,'ScreenSize');
figure('position', [.33*screen([3 4]) .33*screen([3 4])], 'color','w')
% uPo
subplot(121)
title('Prior mean','fontsize',18)
imagesc(uPo)
xlabel('std(likelihood)','fontsize',18)
ylabel('std(prior)','fontsize',18)
set(gca,'fontsize',18)
% sPo
subplot(122)
title('Prior std','fontsize',18)
imagesc(sPo)
xlabel('std(likelihood)','fontsize',18)
ylabel('std(prior)','fontsize',18)
set(gca,'fontsize',18)



%% Another way to visualize the solutions
clear all
screen = get(0,'ScreenSize');
figure('position', [.33*screen([3 4]) .33*screen([3 4])], 'color','w')
ul = 5;
uP = 225;
sPo = 10;
fixeduPo = 115; %wanted output for uPo

pair = 0;
% Loop over all pairs {sp, sl}
for sp = 1 :1: 50
    for sl = 1 :1: 50
        
        % Count iterations
        pair = pair + 1;
        
        % Calculate uPo
        uPo(pair) = (ul/sl^2 + uP/sp^2)*sPo^2;
        
    end 
end

% Draw
hold on
% plot the calculated uPo at each iteration
plot(1 : pair, uPo,'k-',...
    'markersize',3)
% plot the "wanted" output uPo
plot([1 pair],[fixeduPo fixeduPo],'r')
ylim([50 200])
% legend
ylabel('uPo','fontsize',18)
xlabel('pair (sl,sp)','fontsize',18)
set(gca,'fontsize',16)
% legend sPo
text(100,175,strcat('sPo =',num2str(sPo)),'fontsize',18)
% legend "wanted" uPo
text(100,105,strcat('fixed uPo =',num2str(fixeduPo)),'fontsize',18,'color','r')


