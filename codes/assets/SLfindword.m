

% SLfindword.m
%
%     author: steeve laquitaine
%       date: 141123
%    purpose: find words in a sequences of words
%
%usage
%
%           posTarget = SLfindword({'hi','my','friend'},{'hi'})

function posTarget = SLfindword(sentence,targets)

%find position target word in a sentence
if iscell(sentence)    
   
    posTarget = [];
    for i = 1 : length(targets)
        
        posTargetmp = strcmp(sentence,targets{i});
        posTarget = [posTarget; posTargetmp];
        
    end
    
else
    
    %status
    fprintf('%s \n','(SLfindword) Input sentence must be a cell')
    keyboard
end