

%slMake2Dsine
%
%
%A 2D sinusoidal signal is defined
% - its width (w) and height (h), 
% - amplitude (A), 
% - periods (Tw and Th) or frequencies (fw and fh) for width and height 
% - phases (phiw and phih) for width and height.
%
%We can measure the periods in # of pixels (samples), and the
%frequencies in # of cycles per pixel (or cycles per sample). You
%can also specify the periods by n (# of cycles within the signal
%width w and height h).
%
%
%ref:
%http://www.inf.ufsc.br/~visao/khoros/html-dip/c2/s2/front-page.html
%slMake2Dsine(1:100,1:100,1,0.5,0.5,10,10)

function slMake2Dsine(w,h,A,fw,fh,phiw,phih)

s = A*sin(2*pi*fw*w+phiw + 2*pi*fh*h+phih);

imagesc(s)