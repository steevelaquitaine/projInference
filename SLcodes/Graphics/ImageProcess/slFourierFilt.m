
%slFourierFilt.m
%
%
% author: steeve laquitaine
%purpose: fourier filter an image with a low pass mask which means
%         hig frequencies are removed from the image
%         We do low frequency filtering using a Gaussian mask in
%         frequency space.
%
%  date: 151112
% usage:
%
%     im = rand(100,100)>0.5;
%     sigmaf = 10
%     imflt = slFourierFilt(im,sigmaf,'gausWeights')
%
%inputs :
%
%       signaf : gaussian weight filter function standard deviation
%
%Description :
%   1 - The image is fourier transformed
%   2 - Create a gaussian weight filter function
%   3 - multiply the image and the filtered function
%   4 - transform back into the spatial domain
%
%reference :
%http://www.mathworks.com/matlabcentral/fileexchange/28810-fourier-transform-demonstration/content/fourier_demo/html/fourier_demo.html


function imflt = slFourierFilt(im,sigmaf,weights)

%so far we can only keep low frequencies
%make the code so that we can input
%the frequencies we want to keep

% 1 - fourier filtering of the image
tf = fft2(im);

%create reweighting mask
%to filter out high frequencies
[M,N] = size(tf);
mask = zeros(M,N);
[fy,fx] = ndgrid(0:M/2,0:N/2);
ffreq = [M/2 N/2];

%2 - Create a gaussian weight filter function
if strcmp(weights,'gausWeights')

%Create a 2D - Gaussian weight filter function
%The function has its highest values at the four corner of the image
%%keep HFreq (Gaussian mask centred on highest freq)
%mask(1:M/2+1,1:N/2+1) = exp(-((fx-ffreq(2)).^2+(fy-ffreq(1)).^2)/(2*sigmaf)^2);
mask(1:M/2+1,1:N/2+1) = exp(-((fx-1).^2+(fy-1).^2)/(2*sigmaf)^2);%keep low freq
%make symmetrical
mask(1:M/2+1, N:-1:N/2+2) = mask(1:M/2+1, 2:N/2);
mask(M:-1:M/2+2, :) = mask(2:M/2, :);
else
slPrintfStr('slFourierFilt','Input weights... For now weights can only be "gausWeights"')
imflt = [];
return
end

%multiply the image and the filtered function and
%transform back into the spatial domain
imflt = ifft2(mask .* tf);