
% steeve laquitaine 24062009

% save multiple file.fig in the current Directory
% if it says mistake, restart matlab

for i=1:size(action,2)
    plot(action(:,i))
    H=get(0,'children'); %get(root handle,'property');handles of children of the root object (one figure)
    saveas(H,['shyashin' num2str(i) '.fig'])
    clear H
end
