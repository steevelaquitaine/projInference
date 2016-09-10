 %author: steeve laquitaine
   %date: 140117
%purpose: draw half rectified sinusoid
  %usage: channel=drawHalfRectifiedSine(x,peakX)
          %x (e.g., motion direction) and peakX (e.g.,highest response 
          %motion direction) are in radians

function channel=drawHalfRectifiedSine(x,peakX)

%half rectified sine
channel=sin(x+pi/2-peakX);
channel(channel<0)=0;
