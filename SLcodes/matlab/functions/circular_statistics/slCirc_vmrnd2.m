%slCirc_vmrnd2.m
%
% author: steeve laquitaine
%   date: 150906
%purpose: create Nrow by Ncol matrix of values sampled from a von Mises
%         density
%
%  usage:
%
%
%       x = slCirc_vmrnd2(pi,40,250,328);

function x = slCirc_vmrnd2(u,kappa,Nr,Nc)

%svec = -pi:pi/10:pi-pi/10';                            %-pi to pi
svec = SLde2r(5:10:360,0);
v    = 2/pi/besseli(0,kappa,1)*exp(kappa*cos(svec - u) - kappa);   

%case large kappa
if kappa > 1e300
    v   = zeros(1,length(svec));
    v(svec==u) = 1;
end

x1   = randsample(svec,Nc,'true',v);                             %sample v
si   = randi(Nc,Nr,Nc);                                  %sample pos
x    = x1(si);                                                   %Nr by Nc matrix of values
