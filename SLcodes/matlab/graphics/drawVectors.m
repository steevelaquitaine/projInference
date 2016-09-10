%  author: steeve laquitaine
%
%
%    date: 131128
% purpose: draw vectors on a polar
%   usage: 
%   
%           drawVectors([180 90],[1 1]) %for polar angle
%           drawVectors([0 1; 1 0; 0 0]); %for cartesian coordinates
%           the vector's origin is at 0.  
%           "theta" is the angle formed with the vector and the x axis.
%           "length" is the length of the vector.
  
function drawVectors(theta,length)
p=polar(0,1);

%draw vectors
%inputs are angles in degrees
if nargin==2
    for l = 1 : numel(length) 
        coor(l,:) = SLpolar2cartesian(theta(l),length(l),'polar');
    end
    x=coor(:,1);
    y=coor(:,2);
    arrow([0,0],[x y],15,'BaseAngle',60,'EdgeColor',[.5 0.5 0.5],'FaceColor',[.5 0.5 0.5]);    

    %inputs are cartesian coordinates (x,y)
elseif nargin==1
    cartesian=theta;
    arrow([0,0],cartesian,15,'BaseAngle',60,'EdgeColor',[.5 0.5 0.5],'FaceColor',[.5 0.5 0.5]);
end

%get rid of text
%find and remove lines
h=findall(gca,'type','line');
h(h==p)=[];
delete(h);
%find and delete text
t=findall(gca,'type','text');
delete(t);
    