
function Resp=Dchange(Basal,After)

% steeve laquitaine 15042009

% use Basal and After 2 matrix of Basal activity and Response activity for each neurons (columns) and each trials (line)
% output: Resp: Neural response in % change of basal activity


Resp=(After-Basal)
    
