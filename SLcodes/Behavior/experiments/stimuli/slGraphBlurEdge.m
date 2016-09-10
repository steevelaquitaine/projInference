
%slGraphBlurEdge.m
%
% author: steeve laquitaine
%   date: 151216
%purpose: create a matrix of blurred oriented edge area (1) on a black
%         background (0)
%
%usage:
%
%       Wdg = slGraphBlurEdge(400,500,135,100,5)
%
%Description: 
%
%       H: Matrix (image) Height in elements (pix)
%       W: Matrix (image) Width
%     ang: orientation angle in deg (e.g., 135)
%       r: radius in # of pixels
%      sf: sigma of the fourier transform low pass-filter Gaussian mask


function Wdg = slGraphBlurEdge(H,W,ang,r,sf)

%Define the image size and where you want the circle
center = round([W/2 H/2]);  %In [X,Y] pix coord
angMask  = zeros(H,W);

%mark as 1 matrix elements that 
%of the lines that cover all angles
%between ang+10 and ang-10
for i = ang-2:1:ang+2
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

%blur by removing high spatial fr
Wdg = slFourierFilt(angMask,sf,'gausWeights');

%scale between 0 and 1
Wdg = (Wdg - min(Wdg(:)))/(max(Wdg(:)) - min(Wdg(:)));


