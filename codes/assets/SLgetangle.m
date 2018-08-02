
%   authors: steeve laquitaine
%      date: 131128
%     usage: angle=getangle(1,0)
%             theta=arctan(opposite/adjacent);
%             "angle" are in degrees
% reference: http://www.mathsisfun.com/polar-cartesian-coordinates.html           

function angle = SLgetangle(x,y)

%round x and y to 10-4 to prevent weird angle when
%directions cancel each others
x = round(x*10^4)/10^4;
y = round(y*10^4)/10^4;

%convert cartesian to angle (degrees)
angle = atan(y./x); 
deg = SLra2d(angle);

%adjust each angle according to his quadrant (degrees)
for i=1:numel(x)
    
    %(quadrant 1)
    if x(i)>=0 && y(i)>=0              
        angle(i)=angle(i)*180/pi;
    
    %(quadrant 2 & 3)
    elseif x(i)<0                      
        angle(i)=angle(i)*180/pi + 180;
    
    %(quadrant 4)
    elseif x(i)>=0 && y(i)<0           
        angle(i)=angle(i)*180/pi + 360;
    end
end



