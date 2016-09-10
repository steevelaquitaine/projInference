  %Author: Steeve Laquitaine
    %date: 140411
   %usage: matrix=randomizeMatrix(matrix)
 %purpose: randomize a matrix

function matrix=randomizeMatrix(matrix)
matrix(randperm(numel(matrix)))=matrix;
