
% need a matrix of N1, and a matrix of N2
% give a matrix of R1 and a matrix of R2

R1=N1./(N1+N2);
  

for i1= 1:length(N1(1,:))
    for i = 1:length(N1(:,1))
    if N1(i,i1)==0;
        if N2(i,i1)~=0
            R1(i,i1)=0
        end
    end
    if N1(i,i1)~=0
        if N2(i,i1)==0
         R1(i,i1)=1   
        end
    end
    if N1(i,i1)==0
        if N2(i,i1)==0
         R1(i,i1)=0.5   
        end
    end    
    end
end

R2=1-R1;    

