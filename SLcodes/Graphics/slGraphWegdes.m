
%slGraphWegdes.m
%
% author: steeve laquitaine
%   date: 151115
%purpose: create matrix of lines (1) that form a wedge 
%         on a dark background (0)
%
%  usage:
%
%          angMask = slGraphWegdes(135,400,400,100)
%
%Description:
%
%   W: matrix # of rows
%   H: matrix # of cols
%   r: radius wedge
%
%see slGraphWedges2.m



function angMask = slGraphWegdes(ang,W,H,r)

%Define the image size and where you want the circle
center = round([W/2 H/2]);  %In [X,Y] pix coord
angMask  = zeros(H,W);

%mark as 1 matrix elements that 
%of the lines that cover all angles
%between ang+10 and ang-10
for i = ang-10:1:ang+10
    th = SLde2r(i,1);
    xi = round(center(1)+(0:r)*(cos(th)));
    yi = round(center(2)-(0:r)*(sin(th)));
    if max(yi)>H || max(xi)>W
        slPrintfStr('slGraphWegdes',' radius is too large with respect to width and height')
        keyboard
    end
    linIx = sub2ind([H W],yi,xi);
    angMask(linIx) = 1;
end

%draw those lines
%imshow(angMask)