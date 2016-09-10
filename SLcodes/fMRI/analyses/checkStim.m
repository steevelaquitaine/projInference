

%checkStim

cd /Volumes/DroboBKUP/data/mglRetinotopy/s032720160325/Etc

%load data
load('im1.mat')
im1 = im; clear im;
load('im3.mat')
im3 = im; clear im;
load('im5.mat')
im5 = im; clear im;
load('im8.mat')
im8 = im; clear im;

%animate
figure(1)
title('scan 1')
for i = 1 : 477
    subplot(131)
    imagesc(im1(:,:,i)-im3(:,:,i))
    
    subplot(132)
    imagesc(im3(:,:,i))

    subplot(133)
    imagesc(im5(:,:,i))
    
    colormap('gray')
    pause(0.05)
    drawnow
end

%%
figure(2)
subplot(121)
imagesc([im1(:) im3(:)]) 
colormap('gray')

% subplot(122)
% imagesc([im5(:) im8(:)]) 
% colormap('gray')
