
function Window=mobile_wind(bin1,nb_trl)

%%% steeve laquitaine 02/04/2010

%%% create mobile window

%%%  size of mobile window
    %bin1=20;

%%%  size of vector.
    %nb_trl=250;

for i= 1:nb_trl-bin1+1
   Window(i,:)=[i  i+bin1-1];
end

Window(:,1)=[];
