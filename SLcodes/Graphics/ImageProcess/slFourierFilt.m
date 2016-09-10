
%slFourierFilt.m
%
%
%author: steeve laquitaine
%  date: 151112
% usage: 
%       
%       im = rand(10,10)>0.5;      
%       sigmaf = 10
%       imflt = slFourierFilt(im,sigmaf,'gausWeights')


function imflt = slFourierFilt(im,sigmaf,weights)

%so far we can only keep low frequencies
%make the code so that we can input
%the frequencies we want to keep

%fourier filtering of
%low freq
tf = fft2(im);

%create reweighting mask
%to filter out high frequencies
[M,N] = size(tf);
mask = zeros(M,N);
[fy,fx] = ndgrid(0:M/2,0:N/2);
ffreq = [M/2 N/2];

%filter with gaussian weights
if strcmp(weights,'gausWeights')
    
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

%get filtered image
imflt = ifft2(mask .* tf);