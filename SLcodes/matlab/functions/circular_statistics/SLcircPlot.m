%SLcircPlot.m 
%
%  author: Steeve Laquitaine
%    date: 140729
% purpose: draw circular plot
%   usage: figure('color','w'); SLcircPlot([5:10:355],round(vmPdfs([5:10:355],225,2.77,'norm')*1000),5)
% 
% SLcircPlot(1:1:360,vmPdfs(1:1:360,225,2.77,'norm')*1000,5)
% varargin are:
%   - 'scaledRadius'


function SLcircPlot(angles,lengths,baseRadius,iftext,varargin)
hold all

if ~exist('iftext','var')
  iftext=0;
end

%check that angles and lengths are column vectors
if size(angles,1)<size(angles,2)
    angles=angles';
end
if size(lengths,1)<size(lengths,2)
    lengths=lengths';
end

%don't show 0-length arrows
angles(lengths==0)=[];
lengths(lengths==0)=[];

%draw base circle

%calculate x0s & y0s
%We set the radius of the circle that supports the vectors in percent of
%the maximum vector length and draw the base circle.
%scale only if askedso.
if strcmp(varargin,'scaledRadius')
    baseRadius = baseRadius*max(lengths);
end
    
cAx=SLcircle(0,0,baseRadius,[0.75 0.75 0.75]);
axis square
axis off
c0 = SLpolar2cartesian(angles,baseRadius);
x0=c0(:,1);
y0=c0(:,2);

%calculate xEnd & yEnd
c2 = SLpolar2cartesian(angles,lengths);
cEnd=c0+c2;
xEnd=cEnd(:,1);
yEnd=cEnd(:,2);

%draw clean polar & bars
p=polar(0,max(lengths)+baseRadius);
h2=findall(gca,'type','patch');
delete(h2)
h=findall(gca,'type','line');
delete(intersect(h(2:end),p));
t=findall(gca,'type','text');
delete(t);

%radians
rad = SLde2r(angles,1);
polar(rad,lengths+baseRadius);


% arrow([x0 y0],[xEnd yEnd],'Length',0,'BaseAngle',[],'Width',3,'facecolor',...
%     colorH,'edgecolor',colorH)

rlim=max(lengths)+baseRadius;
axis([-1 1 -1 1]*rlim);

%text
if strcmp(iftext,'text')==1
    %lengths
    for i=1:numel(angles)
        if lengths(i)~=0
            text(1.12*xEnd(i),1.12*yEnd(i),num2str(lengths(i)),'fontsize',10,'HorizontalAlignment','center','verticalAlignment','middle')
        end
    end
end
