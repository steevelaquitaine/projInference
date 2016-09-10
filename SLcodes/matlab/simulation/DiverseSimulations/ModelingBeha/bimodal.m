function y3 = bimodal(x,mean1,mean2,stddev1,stddev2)

y1=exp( -((x-mean1).^2)/(2*(stddev1.^2))) / (stddev1 * sqrt(2*pi));
y2=exp( -((x-mean2).^2)/(2*(stddev2.^2))) / (stddev2 * sqrt(2*pi));
y3=0.5*y1+0.5*y2;

y3=y3/sum(y3);
