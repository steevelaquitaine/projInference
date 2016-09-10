%SLdrawhistdir
%
%     author: steeve laquitaine
%       date: 140506
%    purpose: draw custom made histogram
%
function axh = SLdrawhistdir(sample,count,fig,colorH)

%graphics
axh.name = strcat(fig.name,'_his'); 
myBarWidth=1;

%draw
hold all

%histogram
%bar_xend = 2*sample(end)-sample(end-1);
% bar([sample bar_xend],[count count(1)],...
%     'FaceColor',colorH,...
%     'BarWidth',myBarWidth,...
%     'EdgeColor','none')

bar(sample,count,...
    'FaceColor',colorH,...
    'BarWidth',myBarWidth,...
    'EdgeColor','none')

%lines
% plot([sample bar_xend],[count count(1)],...
%     'color',colorH,...
%     'lineWidth',2,...
%     'linesmoothing','on');
plot([225 225],[0 max(count)],'w:')

for j = 1 : numel(sample)
    text(sample(j),count(j),num2str(count(j)),'fontsize',10)
end

%Set x-axis
if numel(sample) > 10
    xunit = 1 : 11 : numel(sample);
else
    xunit = 1 : 1 : numel(sample);
end
set(gca,'fontsize',14,...
    'xtick',sample(xunit),...
    'xticklabel',sample(xunit))
ylabel('Count')
%xlim([min(sample)-10 bar_xend+10])
xlim([0 359])
ylim([0 max(count)]);
box off