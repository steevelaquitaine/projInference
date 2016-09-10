function Sp=drawPrior(fitP,stdPa,Kptrue,coh,s)

figure('color',[1 1 1]);

%Graphics
%subjects colors
colors=[209,213,255;
    189,215,231
    107,174,214
    33,113,181]/255;
for sub=1:s.num
    subplot(1,3,1);
 
    %K
    %..
    %priors' true & fit values
    Sp.es.raw(:,sub)=fitP(4:7,sub);
    if isempty(stdPa)==0;
        Sp.std(:,sub)=stdPa(4:7,sub);
    else
        Sp.std=[];
    end
    
    %draw
    hold all
    myerrorbar(Kptrue,Sp.es.raw(:,sub),'yError',Sp.std(:,sub),...
        'Symbol=o',...
        'Markersize=30',...
        strcat('Color=[',num2str(colors(sub,:)),']'))
    plot(Kptrue,Kptrue,'k:','linewidth',2\3);
    if isnan(Sp.es.raw(:,sub))~=1
        linefit(Kptrue,Sp.es.raw(:,sub)',colors(sub,:))
    end
    title ('Internal prior strength','fontsize',12)
    xlabel('Experimental k_p_r_i_o_r','fontsize',12)
    ylabel('Estimated k_p_r_i_o_r','fontsize',12)
    legend('boxoff')
    ylim([0 1200])
    axis square
    drawPublishAxis
    
    
    %STD (ratio to strong)
    %.....................
    %Prior strength in std (ratio to strong) ;
    Spr(:,sub)=KtoStdev(1:1:360,180',fitP(4:7,sub)',1000);
    Spr(:,sub)=Spr(:,sub)./Spr(4,sub);
    Sper=KtoStdev(1:1:360,180',Kptrue,1000);
    Sper=Sper./Sper(4);
    
    %draw
    subplot(1,3,2);
    hold all
    stdPaSp=[];
    myerrorbar(Sper,Spr(:,sub),'yError',stdPaSp,...
        'Symbol=o',...
        'Markersize=30',...
        strcat('Color=[',num2str(colors(sub,:)),']'));
    if isnan(Spr(:,sub))~=1
        linefit(Sper',Spr(:,sub),colors(sub,:))
    end
    plot(Sper,Sper,'k:','linewidth',3)
    title ({'Internal prior strength','(ratio to strong)'},'fontsize',14)
    xlabel({'Experimental std_p_r_i_o_r','(ratio to strong)'},'fontsize',14)
    ylabel({'Estimated std_p_r_i_o_r', '(ratio to strong)'},'fontsize',14)
    ylim([0 60])
    axis square
    drawPublishAxis
    
    
    %STD
    %...
    %Prior strength in std ;
    Spi(:,sub)=KtoStdev(1:1:360,180',fitP(4:7,sub)',1000);
    Spei=KtoStdev(1:1:360,180',Kptrue,1000);
    
    %draw
    subplot(1,3,3);
    hold all
    stdPaSp=[];
    myerrorbar(Spei,Spi(:,sub),'yError',stdPaSp,...
        'Symbol=o',...
        'Markersize=30',...
        strcat('Color=[',num2str(colors(sub,:)),']'));
    if isnan(Spi(:,sub))~=1
        linefit(Spei',Spi(:,sub),colors(sub,:))
    end
    plot(Spei,Spei,'k:','linewidth',3)
    title ('Internal prior strength','fontsize',14)
    xlabel('Experimental std_p_r_i_o_r','fontsize',14)
    ylabel('Estimated  std_p_r_i_o_r','fontsize',14)
    xlim([0 110])
    ylim([0 110])
    axis square
    drawPublishAxis  
end


% %Get the llh' true and estimated values
% Sl.es.raw=fitP(1:3);
% if isempty(stdPa)==0;
%     Sl.std=stdPa(1:3);
% else
%     Sl.std=zeros(numel(Sl.es.raw),1);
% end
% 
% %Draw subjects' data (i.e.,the strength of the prior as perceived by subjects)
% figure('color',[1 1 1]);
% hold all
% title ('Likelihood strength','fontsize',12)
% xlabel('Coherence','fontsize',12)
% ylabel('Estimated k_l_l_h','fontsize',12)
% myerrorbar(coh,Sl.es.raw,'yError',Sl.std,'Color=[0 0 0]');
% xlim([0 .3])
% ylim([0 max(Sl.es.raw)+max(Sl.std)])
% set(gca,'xtick',0:0.06:.3,'xticklabel',0:0.06:.3,'fontsize',20)
% drawPublishAxis


% %in std instead of k;
% Sli=KtoStdev(1:1:360,180',fitP(1:3)',1000);
% stdPaSl=KtoStdev(1:1:360,180',stdPa(1:3),1000);
% stdPaSl=[0 0 0 0];
% figure('color',[1 1 1]);
% hold all
% myerrorbar(coh,Sli,'yError',stdPaSl,...
%     'Symbol=o',...
%     'Markersize=30',...
%     'Color=[0 0 0]');
% title ('Likelihood strength','fontsize',12)
% xlabel('Coherence','fontsize',12)
% ylabel('Estimated std_l_l_h','fontsize',12)
% %set(gca,'xtick',0:10:90,'xticklabel',0:10:90,'fontsize',20)
% drawPublishAxis
