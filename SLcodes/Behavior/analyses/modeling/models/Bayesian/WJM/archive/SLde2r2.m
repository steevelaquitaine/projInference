%SLde2r2.m
%
%       $Id: SLde2r2.m 750 2012-09-08 07:12:46Z steeve $
%     usage: SLde2r2(180,1)
%        by: steeve laquitaine
%      date: 22/10/12 last modification: 140610 21:34:00


function radians = SLde2r2(ang,sign)

%call for help (may slow down code)
% if ieNotDefined('ang')
%     help de2r
%     return
% end

%memory and speed
ang = single(ang);

%not signed radians (1:2*pi)
radians = (ang/360)*2*pi;

%sign radians(-pi:pi)
if sign==1
    radians(ang>180)=(ang(ang>180)-360)*(2*pi/360);
end