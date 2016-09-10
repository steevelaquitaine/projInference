%SLcalculateVectNorm.m
%
%       $Id: SLcalculateVectNorm.m 750 2012-09-08 07:12:46Z steeve $
%     usage: Vectnorm=SLcalculateVectNorm(x,y)
%        by: steeve laquitaine
%      date: 140526
%   purpose: calculate the norm of a vector: pythagoras formula

function Vectnorm=SLcalculateVectNorm(x,y)

Vectnorm=sqrt(x^2+y^2);

