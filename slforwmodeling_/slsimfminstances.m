
%slsimfminstances.m
%
%
% author : steeve laquitaine
%purpose : generates S response instances x N voxels
%          for a vector of S stimulus directions
%          stimulus directions values must belong to 1:1:360 deg
%
%usage :
%
%         svec = [15*ones(10,1); 85*ones(10,1)];
%         W = rand(100,8); 
%         [b_wres,b,fm] = slsimfminstances(svec,W,sigma,tau,rho)
%
%inputs
%   
%          svec : vector of trial instances of stimulus directions
%             W : a matrix of voxels x channel positive weights
%
%outputs:
%
%             
%        b_wres : a matrix of S response instances x N voxels %                 
%             C : a matrix of the channel average responses to the
%                 vector of stimulus directions
%                 Ni instances x 8 channels
%             b : a matrix of S response instances x N voxels 
%            fm : structure containing model architecture and 
%                 parameters
%

function [b,C,fm] = slsimfminstances(svec,W)


%5 direction tuned channels
fm = slsimfmchannels;
f_k_s = fm.f_k_s';

%average channel responses to vector of stimulus directions
%Ni instances x 8 channels
C = fm.f_k_s(:,svec)';

%simulate a matrix b of Ni trial stimulus instances x N voxels 
Ni = length(svec);
N = size(W,1);
b = nan(Ni,N);
for i = 1 :  N   
    Wi = W(i,:);
    b(:,i) = sum(Wi(ones(1,Ni),:).*f_k_s(svec,:),2);
end

