
%SLMixtureVM.m
%
% author: Steeve Laquitaine
%   date: 150109
%purpose: Create a mixture of two von Mises proability distribution
%
%
%usage:
%
%       task = SLmakeDiscreteMixtureVM(8.4,5:10:355,[145 305],200)
%  

%mixture of Von Mises (MoVM)
function MoVM = SLMixtureVM( x, modes, k)

MoVM = 0.5*vmPdfs(x,modes(1),k,'norm') + 0.5*vmPdfs(x,modes(2),k,'norm');
