function [ Sig,PValue,ContigenMatrix ] = FisherExactTest( XVector,YVector )

% This function takes into account of n*m contingency table, instead of 2*2
% or 3*3.
% The idea of Fisher's Exact Test here is learned from website:
% http://mathworld.wolfram.com/FishersExactTest.html

% If Sig=1, then variables in XVector and YVector are significantly associate.
% If Sig=0, then variables in XVector and YVector are not associate.

% Assumption for sample XVector and YVector
% 1. both of them have the same length
% 2. missing data is represented as "-1", of course, one can change it in variable "Miss"
% 3. assume the confidence of test is 0.05, one can change it in "Conf"

% This function was written by Lowell Guangdi at 2009/06/08,revised at
% 2010/01/28

Conf = 0.05; % significant level is set to 0.05, users can change it.
Miss = -1;
if nargin > 1  % Input two vectors representing two variables.
% If two variables do not have the same input length of sample, give a
% warning
if length( YVector ) ~= length( YVector )
   PValue = NAN;  
   return
end
N  = length( YVector );
XState = unique( XVector );
YState = unique( YVector );

% Select out the distinguished state
if isequal( XState( 1 ), Miss ) == 1  
    XState( 1 )= [];
end
if isequal( YState( 1 ), Miss ) == 1
    YState( 1 )= [];
end

% Find out the contingency table
ContigenMatrix = zeros( length( XState), length( YState ) );
for p = 1:N
    if isequal( XVector( p ),Miss ) == 0 && isequal( YVector( p ),Miss ) == 0 
        Row = find( XState == XVector( p ),1);
        Col = find( YState == YVector( p ),1);
        ContigenMatrix( Row,Col ) = ContigenMatrix( Row,Col ) + 1;
    end
end
else
    ContigenMatrix = XVector;
end

% With contingency table , compute the p-value according the definition 
TempTable = ContigenMatrix;
PValue = 1;
RowVec = sum( TempTable,2 ); ColVec = sum( TempTable,1 );
LenRow = length( RowVec ); LenCol = length( ColVec );
n = sum( sum( TempTable ) );
Run1 = 1; Run2 = 1; Count = 0;
while ~( Run1 == 0 && Run2 == 0 )
   Count = Count + 1;
   % Consider to divide the elements in aij!
   p = 1; Flag = 1;
   while Run1 == 1 && Flag == 1 && p <= LenRow
       q = 1;
       while Flag == 1 && q <= LenCol
             if PValue < 0.00001
                Flag = 0;
             elseif TempTable( p,q ) > 1 
                PValue = PValue / TempTable( p,q );
                TempTable( p,q ) = TempTable( p,q )-1;
             elseif TempTable( p,q ) <= 1 
                 q = q + 1;               
             end            
       end
       p = p + 1;
   end
   % consider to divide the elements in n!
   while Run1 == 1       
       if n == 1 || PValue < 0.00001
          break;
       else  %if n > 1 && PValue > 0.00001
          PValue = PValue / n;
          n = n - 1;
       end
   end
   % if all the elements in n! and aij are completely divided.
   if n==1 && p > LenRow && q > LenCol
       Run1 = 0;
   end
   % consider to multiply the elements in Ri and Ci
   p = 1; q = 1; 
   while Run2 == 1
       % consider to multiply Ri
       if PValue > 10000
           break;
       elseif p <= LenRow && RowVec( p ) > 1 
          PValue = PValue * RowVec( p ); 
          RowVec( p ) = RowVec( p ) - 1;
       elseif p <= LenRow %&& RowVec( p ) <= 1
           p = p + 1;
       end
       % consider to divide by Ci
       if PValue > 10000
           break;
       elseif q <= LenCol && ColVec( q ) > 1
          PValue = PValue * ColVec( q ); 
          ColVec( q ) = ColVec( q ) - 1;
       elseif q <= LenCol %&& ColVec( q ) <= 1
           q = q + 1;
       else
       end
       % if all the elements of Ci! and Ri! are divided.
       if  q > LenCol && p > LenRow
           Run2 = 0;
       end       
   end
   % For sure that the p-value will be less than given confidence.
   if PValue < Conf && Run1 == 0 || Count > 100000
       break;
   end   
end
% assume the significant confidence is 0.05, if p-value is less than 0.05
% ,then it is safely to say the variable in XVector and YVector are
% significant associated --that is, they are highly dependent. 
Sig = 0;
if PValue <= Conf  
    Sig = 1;
end
end

