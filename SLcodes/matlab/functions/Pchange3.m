
function Resp=Pchange2(Basal,After)

% steeve laquitaine 15042009

% use Basal and After 2 matrix of Basal activity and Response activity for each neurons (columns) and each trials (line)
% output: Resp: Neural response in % change of basal activity



% core of Pchange
for i=1: length(Basal(:,1))
    EBasal(i,:)=nanmean(Basal)
end
Resp=(After-Basal)./Basal;

    
% for Resp=NaN (0/0 or 0/nb))
for i2=1:length(Basal(1,:))
    for i3=1:length(Basal(:,1))
        if After(i3,i2)-Basal(i3,i2)==0
            Resp(i3,i2)=0;
        end
    end
end

% for Resp=infinite (nb/0)
for i2=1:length(Basal(1,:))
    for i3=1:length(Basal(:,1))
        if Basal(i3,i2)==0
        Resp(i3,i2)=After(i3,i2)-Basal(i3,i2)
        end
    end
end

