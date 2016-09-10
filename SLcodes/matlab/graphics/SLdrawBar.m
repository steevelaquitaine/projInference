    
%SLdrawBar.m
%
%      author: steeve laquitaine
%        date: 140207 last modification: 140513
%     purpose: plot bars or errorbars
%
%       usage: 
%
%           SLdrawBar([1 2 3 4],{'steeve','John','Celine','Minori'},[1 2 4 5])
%
% 
% varargins:
%
%                   'yError', error :
%            'ColorShading',shading : shading is scalar between 0 and 1
%              'FaceColor',facecolor: e.g., facecolor is [.5 .5 .5]

function SLdrawBar(y,names,x,varargin)

%get shading
if any(strcmp(varargin,'ColorShading'))
    shading=varargin{find(strcmp(varargin,'ColorShading'))+1};
    colorf=repmat((0:1/length(y):0.95)'.^shading,1,3);
else
    shading=0.5;
    colorf = repmat((0:1/length(y):0.95)'.^shading,1,3);
end

%or just color
if any(strcmp(varargin,'FaceColor'))
    colorf = repmat(varargin{find(strcmp(varargin,'FaceColor'))+1},length(y),1);
end

%just bar
%--------
for i=1:numel(x) 
    hold all
    b(i)=bar(x(i),y(i),...
        'edgecolor','none',...
        'facecolor',colorf(i,:));
    
    %annotate
    if ~sum(strcmp(varargin,'text'))
        
        text(x(i),1.1*y(i),num2str(fix(y(i)*100)/100),...
            'HorizontalAlignment','center','fontsize',11)
    end
    
    %case text is input
    if sum(strcmp(varargin,'text'))
        textpos = find(strcmp(varargin,'text'))+1;
        texti   = varargin{textpos};
        %annotate
        text(x(i),1.1*y(i),num2str(texti(i)),...
            'HorizontalAlignment','center','fontsize',11)
    end
end

%case errorbar
%-------------
if sum(strcmp(varargin,'yError'))
    
    %get error 
    error=varargin{find(strcmp(varargin,'yError'))+1};
    
    %plot errorbar
    for i=1:numel(x) 
        SLerrorbar(x(i),y(i),'yError',error(i),'MarkerSize=1',...
           'FaceColor',colorf(i,:))
    end
end

%Graphics details
box off
set(gcf,'color','w')

%case names are input
if ~isempty(names)
    set(gca,'fontsize',14,'xtick',x(1:numel(y)),'xticklabel',names)
    rotateXLabels(gca(),45)
else
    set(gca,'fontsize',14,'xtick',x(1:numel(y)),'xticklabel',[])
end






