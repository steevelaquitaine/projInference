% The Mack-Skillings Statistical Test
% A nonparametric two-way ANOVA used for unbalanced incomplete block 
% designs, when the number of observations in each treatment/block pair is
% one or greater.  The test is equivalent to the Friedman test when
% balanced and there are no missing observations.
% 
% Syntax:
% [p stats] = mackskill(response, treatments, blocks)
% [p stats] = mackskill(M, reps)
%   This version assumes M is structured similar to a table in the function
%   friedman(), where columns are treatments (k) (the variable of interest), 
%   and rows are blocks (n) (a nuisance variable).  Reps (c)
%   corresponds to the maximum number of repeated observations.  
%   The M matrix should therefore have k columns and n*c rows.
%
% Example:
%
% The friedman() function uses popcorn data originally published by Hogg
% (1987).  From the text of that function:
% "The columns of the matrix popcorn are brands (Gourmet, National,
% and Generic). The rows are popper type (Oil and Air). The study popped 
% a batch of each brand three times with each popper. The values are the 
% yield in cups of popped popcorn."  Given that the data is structured this
% way, we must be asking "Which brand of popcorn pops the best, apart from
% the type of popper it is popped in."  k=3,n=2,c=3.
%
% load popcorn
% [p table stats] = friedman(popcorn,3);
% % chisq = 13.76
% % p = 0.001
% 
% If there are missing observations, use the Mack-Skilling test:
%
% popcorn(5) = NaN;
% [p stats] = mackskill(popcorn,3);
% stats.T = 8.5856
% stats.p = 0.0137
%
% 
% Example:
%
% Hollander and Wolfe (1999) describe the Mack-Skillings test in their
% example 7.9 "Determination of Niacin in Bran Flakes" on page 331, from
% data originally published in Campbell & Pelletier (1962).
%
% % Table 7.20  Amount of Niacin in Enriched Bran Flakes 
% %
% %         Enrichment (mg/100g of Bran Flakes)
% %          1            4           8
% M =       [7.58       11.63       15.00; ...  % lab 1
%           7.87        11.87       15.92; ...
%           7.71        11.40       15.58; ...
%           8.00        12.20       16.60; ... % lab 2
%           8.27        11.70       16.40; ...
%           8.00        11.80       15.90; ...
%           7.60        11.04       15.87; ... % lab 3
%           7.30        11.50       15.91; ...
%           7.82        11.49       16.28; ...
%           8.03        11.50       15.10; ... % lab 4
%           7.35        10.10       14.80; ...
%           7.66        11.70       15.70];
%
%
% However, we want to know whether the _labs_ are different, apart from the
% amount of niacin enrichment, or as Hollander and Wolfe put it, "Of
% primary interest here is the precision of the laboratory procedure for
% determining niacin content in bran flakes."  Our table is therefore not
% formatted correctly to call mackskill(M,3), as this would tell us whether
% our columns were notably different, apart from rows.
%
% We can reformat our data this way:
%
% [row enrichment niacin] = find(sparse(M)); 
% lab = ceil(row/3);  % 3, since there are 3 repetitions
% response = niacin; treatment = lab; block = enrichment;
%
% [p stats] = mackskill(response,treatment,block)
%
% stats.T  = 12.9274
% stats.df = 3
% stats.p  = 0.0048
%
% These results are equivalent to the test statistic in Hollander & Wolfe
% (1999), MS = 12.93 (page 332).
%
%
% The algorithm uses the chi-squared approximation for the p-value, which 
% should not be used when there are very few observations.  Please refer to
% the original text for a complete description.
% 
% References: 
% Hollander, M., & Wolfe, D. A. (1999). Nonparametric statistical methods (2nd ed.). New York: Wiley.
% Mack, G. A., & Skillings, J. H. (1980). A Friedman-Type Rank Test for Main Effects in a 2-Factor Anova. Journal of the American Statistical Association, 75(372), 947-951.
% Skillings, J. H., & Mack, G. A. (1981). On the Use of a Friedman-Type Statistic in Balanced and Unbalanced Block Designs. Technometrics, 23(2), 171-177.
%
% The code was tested against several published datasets.  For a copy of
% these datasets and other code written by the author, please see:
% http://www.geog.ucsb.edu/~pingel/matlabCode/index.html
%
% Use of this code for any non-commercial purpose is granted under the GNU
% Public License.  
%
% Author: 
% Thomas J. Pingel
% Department of Geography
% University of California, Santa Barbara
% 11 November 2010

function [p stats rankedobs] = mackskill(M, treatments,blocks)
% load popcorn
% M = popcorn;
% reps = 3;

if nargin<3
% This section reformats matrix M into:
%     X (observations in the matrix M)
%     treatments (columns of M, and the variable of interest) 
%     blocks (rows of M, and the nuisance variable)
%     Second argument is assumed to be reps

% Pick apart the observations
reps = treatments;
X = reshape(M,numel(M),1); 

% Since input is a matrix, define the levels.
treatmentlevels = [1:size(M,2)]; % Columns
blocklevels = size(M,1)/reps; % Rows
blocklevels = [1:blocklevels];
% treatmentlevels = ([1:size(M,2)]'); % Columns
% blocklevels = ([1:size(M,1)]'); % Rows
k = length(treatmentlevels);
n = length(blocklevels);
% Create a vector of treatment and observation levels.
treatmentsMatrix = repmat(treatmentlevels,n*reps,1);
treatments = reshape(treatmentsMatrix,numel(X),1); 
blocksMatrix = repmat(reshape(repmat(blocklevels,reps,1),n*reps,1),1,k);
blocks = reshape(blocksMatrix,numel(X),1);
% blocks = reshape(repmat(reshape(repmat([1:n],k,1),n*reps,1),1,k),numel(M),1);

stats.treatmentsMatrix = treatmentsMatrix;
stats.blocksMatrix = blocksMatrix;

%% For testing, keep these handy in double form
tr = treatments;
br = blocks;
% Get rid of extraneous information, as this will be redefined in the next
% section anyway.
clear treatmentlevels blocklevels k n;
end
%%
% This section applies to if the preferred format is supplied
% skillmack(X,treatments,blocks) where X is a vector (double) and
% treatments and blocks are cell arrays.

% X is now the first argument, input as M

if nargin==3
    X = M;
end

% First, convert to a cell array from matrix if necessary
if ~iscell(treatments)
    treatments2 = cell(size(treatments));
    for i=1:length(treatments)
        treatments2{i,1} = treatments(i);
    end
    treatments = treatments2;
    clear treatments2 i;
end
if ~iscell(blocks)
    blocks2 = cell(size(blocks));
    for i=1:length(blocks)
        blocks2{i,1} = blocks(i);
    end
    blocks = blocks2;
    clear blocks2 i;
end
%%
% Change to cell array of strings, for standardization.
for i=1:length(blocks)
    blocks{i,1} = num2str(blocks{i});
    treatments{i,1} = num2str(treatments{i});
end  
clear i;

%% Remove NaNs from cell array

indx = find(isnan(X));
X(indx) = [];
blocks(indx) = [];
treatments(indx) = [];



%%
% Determine unique levels
treatmentlevels = unique(treatments);
blocklevels = unique(blocks);


%%
% Check to see if any block has only one observation.  If so, for now just
% issue a warning.  Technically, this block should be removed.
% [this isn't the case for mack-skill, but leave in for later
% errorchecking]
% for i=1:length(blocklevels)
%     indx = find(strcmp(blocks,blocklevels{i}));
%     if length(indx) <= 1
%         disp(['Block ',num2str(blocklevels{i}),' has an insufficient number of observations.']);
%     end
% end
% clear i indx;

%% Balance the observations
% See if the results improve if 'unbalanced' setups are replaced with NaNs


%%
% Create a vector to hold ranked observations
rankedobs = nan(size(X));
% disp(num2str(size(X)));
% disp(num2str(length(blocklevels)));
for i=1:length(blocklevels)
   % Step II
   % Within each block, rank the observations from 1 to ki, where ki is the
   % number of treatments present in block i.  If ties occur, use average
   % ranks.
   % Grab the blocks at level i
   indx = find(strcmp(blocks,blocklevels{i}));
   % r holds the ranks for that block. NaNs in empty values.
   r = tiedrank(X(indx));
   
   for j=1:length(indx)
      rankedobs(indx(j)) = r(j);
   end
end
clear i j indx indx2 r replacementr
% disp(num2str(rankedobs));
stats.rankedobs = rankedobs;
% stats.rankedobs2 = reshape(rankedobs,size(M));
%% Sum the ranks for each treatment and block
I = length(blocklevels);
J = length(treatmentlevels);
rankblock = nan(I,J);
for i=1:I % rows, blocks
    for j=1:J % treatments, columns
        indx = intersect(find(strcmp(treatments,treatmentlevels{j})),find(strcmp(blocks,blocklevels{i})));
        rankblock(i,j) = sum(rankedobs(indx));
    end
end
clear i j indx;
stats.rankblock = rankblock;
% stats.rankblock2 = reshape(rankblock,size(M));
%% Modify these by treatment
I = length(blocklevels);
J = length(treatmentlevels);
Rj = nan(J,1);
for j=1:J  % column, treatment
   s = 0;
   for i=1:I  % row, blocks
      s = s + (rankblock(i,j) / length(find(strcmp(blocks,blocklevels{i}))));
   end
   Rj(j) =  s;
   clear s;
end
R = Rj - mean(Rj);
clear s j i
        

%% Calculate sigma
I = length(blocklevels);
J = length(treatmentlevels);
sigma = nan(J,J);
for t=1:J  % treatments, columns
    for j=1:J  % also treatments, columns
        s = 0;
        for i = 1:I  % blocks, rows
            nij = length(intersect(find(strcmp(treatments,treatmentlevels{j})),find(strcmp(blocks,blocklevels{i}))));
            nit = length(intersect(find(strcmp(treatments,treatmentlevels{t})),find(strcmp(blocks,blocklevels{i}))));
            ni = length(find(strcmp(blocks,blocklevels{i})));
            
            s = s + ((nij*nit*(ni+1))/(12*(ni^2)));
%             clear nij nit ni;
        end
        sigma(t,j) = -s;
%         clear s
    end
end
clear t j i

for i=1:J
    sigma(i,i) = 0;
    sigma(i,i) = abs(sum(sigma(i,:)));
end

%% Let's try step 4: Calculating weights.
% A = nan(length(treatmentlevels),1);
% maxrank = nan(size(rankedobs));
% frontweight = nan(size(rankedobs));
% backweight = nan(size(rankedobs));
% totalweight = nan(size(rankedobs));
% % Calculate front and back weights 
% for i=1:numel(X)
%    maxrank(i) = max(rankedobs(find(strcmp(blocks,blocks{i}))));
%    frontweight(i) = sqrt(12/(maxrank(i)+1));
%    backweight(i) = rankedobs(i) - ((maxrank(i) + 1)/2);
% end
% 
% % Multiply them together to get total weights
% totalweight = frontweight.*backweight;
% % Sum each treatment.
% for i=1:length(A)
%     indx = find(strcmp(treatments,treatmentlevels{i}));
%     A(i) = sum(totalweight(indx));
% end
% clear i totalweight frontweight backweight maxrank indx;
% % disp(num2str(A));
% 
% %% Create sigma matrix
% sigma = nan(length(treatmentlevels),length(treatmentlevels));
% k = length(treatmentlevels);
% for i=1:k % row
%     for j=1:k % column
%        indxi = intersect(find(strcmp(treatments,treatmentlevels{i})),find(isfinite(X)==1));
%        indxj = intersect(find(strcmp(treatments,treatmentlevels{j})),find(isfinite(X)==1));
% %        indxk = intersect(indxi,indxj);
%        sigma(i,j) = -length(intersect([blocks{indxi}],[blocks{indxj}]));
%     end
% end
% for i=1:length(treatmentlevels)
%     j = setdiff([1:length(treatmentlevels)],i);
%     sigma(i,i) = sum(abs(sigma(i,j)));
% end
% 

%% Calculate the final statistic.
T = R' * pinv(sigma) * R;
df = length(treatmentlevels)-1;
p = 1 - chi2cdf(T,df);
stats.T = T;
stats.df = df;
stats.p = p;
stats.labels = treatmentlevels;
stats.meanranks = Rj';
stats.source = 'mackskill';
stats.n = df;
stats.rankblock = rankblock;

stats.X = X;
stats.blocks = blocks;
stats.treatments = treatments;
