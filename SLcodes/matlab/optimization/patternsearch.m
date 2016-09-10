
%ref: https://lost-contact.mit.edu/afs/cs.stanford.edu/package/matlab-r2009b/matlab/r2009b/toolbox/gads/gads/patternsearch.m

function [X,FVAL,EXITFLAG,OUTPUT] = patternsearch(FUN,initialX,Aineq,Bineq,Aeq,Beq,LB,UB,nonlcon,options)
%PATTERNSEARCH Constrained optimization using pattern search.
%   PATTERNSEARCH attempts to solve problems of the form:
%       min F(X)  subject to:  A*X  <= B, Aeq*X  = Beq (linear constraints)
%        X                     C(X) <= 0, Ceq(X) = 0 (nonlinear constraints)
%                              LB <= X <= UB
%
%   X = PATTERNSEARCH(FUN,X0) starts at X0 and finds a local minimum X to
%   the function FUN. FUN accepts input X and returns a scalar function
%   value evaluated at X. X0 may be a scalar or vector.
%
%   X = PATTERNSEARCH(FUN,X0,A,b) starts at X0 and finds a local minimum X
%   to the function FUN, subject to the linear inequalities A*X <= B.
%
%   X = PATTERNSEARCH(FUN,X0,A,b,Aeq,beq) starts at X0 and finds a local
%   minimum X to the  function FUN, subject to the linear equalities
%   Aeq*X = Beq as well as A*X <= B. (Set A=[] and B=[] if no inequalities
%   exist.)
%
%   X = PATTERNSEARCH(FUN,X0,A,b,Aeq,beq,LB,UB) defines a set of lower and
%   upper bounds on the design variables, X, so that a solution is found in
%   the range LB <= X <= UB. Use empty matrices for LB and UB if no bounds
%   exist. Set LB(i) = -Inf if X(i) is unbounded below;  set UB(i) = Inf if
%   X(i) is unbounded above.
%
%   X = PATTERNSEARCH(FUN,X0,A,b,Aeq,beq,LB,UB,NONLCON) subjects the
%   minimization to the constraints defined in NONLCON. The function
%   NONLCON accepts X and returns the vectors C and Ceq, representing the
%   nonlinear inequalities and equalities respectively. PATTERNSEARCH
%   minimizes FUN such that C(X)<=0 and Ceq(X)=0. (Set LB=[] and/or UB=[]
%   if no bounds exist.)
%
%   X = PATTERNSEARCH(FUN,X0,A,b,Aeq,beq,LB,UB,NONLCON,options) minimizes
%   with the default optimization parameters replaced by values in the
%   structure OPTIONS. OPTIONS can be created with the PSOPTIMSET function.
%   See PSOPTIMSET for details.
%
%   X = PATTERNSEARCH(PROBLEM) finds the minimum for PROBLEM. PROBLEM is a
%   structure that has the following fields:
%      objective: <Objective function>
%             x0: <Starting point>
%          Aineq: <A matrix for inequality constraints>
%          bineq: <B vector for inequality constraints>
%            Aeq: <A matrix for equality constraints>
%            beq: <B vector for equality constraints>
%             lb: <Lower bound on X>
%             ub: <Upper bound on X>
%        nonlcon: <Nonlinear constraint function>
%        options: <options structure created with PSOPTIMSET>
%         solver: <solver name 'patternsearch'>
%       rngstate: <State of the random number generator>
%   This syntax is specially useful if you export a problem from
%   PSEARCHTOOL and use it from the command line to call PATTERNSEARCH.
%   NOTE: PROBLEM must have all the fields as specified above.
%
%   [X,FVAL] = PATTERNSEARCH(FUN,X0,...) returns FVAL, the value of the
%   objective function FUN at the solution X.
%
%   [X,FVAL,EXITFLAG] = PATTERNSEARCH(FUN,X0,...) returns EXITFLAG which
%   describes the exit condition of PATTERNSEARCH. Possible values of
%   EXITFLAG and the corresponding exit conditions are
%
%     1  Magnitude of mesh size is less than specified tolerance and
%         constraint violation less than options.TolCon.
%     2  Change in X less than the specified tolerance and
%         constraint violation less than options.TolCon.
%     3  Change in FVAL less than the specified tolerance and
%         constraint violation less than options.TolCon.
%     4  Magnitude of step smaller than machine precision and
%         constraint violation less than options.TolCon.
%     0  Maximum number of function evaluations or iterations reached.
%    -1  Optimization terminated by the output or plot function.
%    -2  No feasible point found.
%
%   [X,FVAL,EXITFLAG,OUTPUT] = PATTERNSEARCH(FUN,X0,...) returns a
%   structure OUTPUT with the following information:
%          function: <Objective function>
%       problemtype: <Type of problem> (Unconstrained, Bound constrained or
%                     linear constrained)
%        pollmethod: <Polling technique>
%      searchmethod: <Search technique used>, if any
%        iterations: <Total iterations>
%         funccount: <Total function evaluations>
%          meshsize: <Mesh size at X>
%     maxconstraint: <Maximum constraint violation>, if any
%           message: <PATTERNSEARCH termination message>
%
%   Examples:
%    FUN can be a function handle (using @)
%      X = patternsearch(@lincontest6, ...)
%    In this case, F = lincontest6(X) returns the scalar function
%    value F of the  function  evaluated at X.
%
%   An example with inequality constraints and lower bounds
%    A = [1 1; -1 2; 2 1];  b = [2; 2; 3];  lb = zeros(2,1);
%    [X,FVAL,EXITFLAG] = patternsearch(@lincontest6,[0 0],A,b,[],[],lb);
%
%     FUN can also be an anonymous function:
%        X = patternsearch(@(x) 3*sin(x(1))+exp(x(2)),[1;1],[],[],[],[],[0 0])
%     returns X = [0;0].
%
%   If FUN or NONLCON are parameterized, you can use anonymous functions to
%   capture the problem-dependent parameters. Suppose you want to minimize
%   the objective given in the function myobj, subject to the nonlinear
%   constraint myconstr, where these two functions are parameterized by
%   their second argument a1 and a2, respectively. Here myfit and myconstr
%   are M-file functions such as
%
%        function f = myobj(x,a1)
%        f = exp(x(1))*(4*x(1)^2 + 2*x(2)^2 + 4*x(1)*x(2) + 2*x(2) + a1);
%
%   and
%
%        function [c,ceq] = myconstr(x,a2)
%        c = [1.5 + x(1)*x(2) - x(1) - x(2);
%              -x(1)*x(2) - a2];
%        % No nonlinear equality constraints:
%         ceq = [];
%
%   To optimize for specific values of a1 and a2, first assign the values
%   to these two parameters. Then create two one-argument anonymous
%   functions that capture the values of a1 and a2, and call myobj and
%   myconstr with two arguments. Finally, pass these anonymous functions to
%   PATTERNSEARCH:
%
%     a1 = 1; a2 = 10; % define parameters first
%     options = psoptimset('Display','iter'); % Display iterative output
%     x = patternsearch(@(x)myobj(x,a1),[1;2],[],[],[],[],[],[],@(x)myconstr(x,a2),options)
%
%   See also PSOPTIMSET, GA, PSOUTPUTFCNTEMPLATE, SEARCHFCNTEMPLATE, @.

%   Copyright 2003-2008 The MathWorks, Inc.
%   $Revision: 1.20.4.18 $  $Date: 2008/06/20 07:56:11 $

% Old syntax X = PATTERNSEARCH(FUN,X0,A,b,Aeq,beq,LB,UB,OPTIONS) may work
% fine but users are encouraged to update to new syntax of PATTERNSEARCH
% which takes ninth argument as NONLCON and OPTIONS is passed as tenth
% argument.
defaultopt = struct('TolMesh', 1e-6, ...
    'TolCon', 1e-6, ...
    'TolX', 1e-6 , ...
    'TolFun',1e-6 , ...
    'TolBind',1e-3, ...
    'MaxIter', '100*numberofvariables', ...
    'MaxFunEvals', '2000*numberofvariables', ...
    'TimeLimit', Inf, ...
    'MeshContraction', 0.5, ...
    'MeshExpansion', 2.0, ...
    'MeshAccelerator','off', ...
    'MeshRotate','on', ...
    'InitialMeshSize', 1.0, ...
    'ScaleMesh', 'on', ...
    'MaxMeshSize', inf, ...
    'InitialPenalty', 10, ...
    'PenaltyFactor', 100, ...
    'PollMethod', 'gpspositivebasis2n', ...
    'CompletePoll','off', ...
    'PollingOrder', 'consecutive', ...
    'SearchMethod', [], ...
    'CompleteSearch','off', ...
    'Display', 'final', ...
    'OutputFcns', [], ...
    'PlotFcns',[], ...
    'PlotInterval', 1, ...
    'Cache', 'off', ...
    'CacheSize',1e4, ...
    'CacheTol',eps, ...
    'Vectorized','off', ...
    'UseParallel','never' ...
    );

% If just 'defaults' passed in, return the default options in X
if nargin == 1 && nargout <= 1 && isequal(FUN,'defaults')
    X = defaultopt;
    return
end

errmsg = nargchk(1,10,nargin);
% At least 1 arguments are needed.
if ~isempty(errmsg)
    error('gads:patternsearch:inputArg',[errmsg,' PATTERNSEARCH requires at least 1 input argument.']);
end

if nargin < 10,  options = [];
    if nargin < 9,  nonlcon = [];
        if nargin < 8, UB = [];
            if nargin < 7, LB = [];
                if nargin <6, Beq = [];
                    if nargin <5, Aeq = [];
                        if nargin < 4, Bineq = [];
                            if nargin <3, Aineq= [];
                            end
                        end
                    end
                end
            end
        end
    end
end

% Ninth argument (nonlcon) could be options structure (old syntax from ver 1.0)
if isstruct(nonlcon) && nargin < 10
    options = nonlcon;
    nonlcon = [];
end

% One input argument is for problem structure
if nargin == 1
    if isa(FUN,'struct')
        [FUN,initialX,Aineq,Bineq,Aeq,Beq,LB,UB,nonlcon,rngstate,options] = separateOptimStruct(FUN);
        % Reset the random number generators
        resetDfltRng(rngstate);
    else % Single input and non-structure.
        error('gads:patternsearch:inputArg','The input should be a structure with valid fields or provide at least two arguments to PATTERNSEARCH.' );
    end
end

try
    dataType = superiorfloat(initialX,Aineq,Bineq,Aeq,Beq,LB,UB);
    if ~isequal('double', dataType)
        error('gads:patternsearch:dataType', ...
            'PATTERNSEARCH only accepts inputs of data type double.')
    end
catch
    error('gads:patternsearch:dataType', ...
        'PATTERNSEARCH only accepts inputs of data type double.')
end

% If FUN is a cell array with additional arguments, handle them
if iscell(FUN)
    objFcnArg = FUN(2:end);
    FUN = FUN{1};
else
    objFcnArg = {};
end

% Only function_handle or inlines are allowed
if isempty(FUN) ||  ~(isa(FUN,'inline') || isa(FUN,'function_handle'))
    error('gads:patternsearch:needFunctionHandle','Objective function must be a function handle.');
end

% If NONLCON is a cell array with additional arguments, handle them
if iscell(nonlcon)
    conFcnArg = nonlcon(2:end);
    nonlcon = nonlcon{1};
else
    conFcnArg = {};
end
% Constraint function must be a function_handle or inline
if ~isempty(nonlcon) && ~(isa(nonlcon,'inline') || isa(nonlcon,'function_handle'))
    error('gads:patternsearch:needFunctionHandle','Constraint function must be a function handle.');
end

if(~isempty(initialX))
    X = initialX;
    Iterate.x = initialX(:);
    numberOfVariables = length(Iterate.x);
else
    error('gads:patternsearch:initialPoint','You must provide an initial point.');
end

% Determine the 'type' of the problem
if ~isempty(nonlcon)
    type = 'nonlinearconstr';
    % Determine the sub-problem type for the constrained problem (used in ALPS)
    if ~isempty(Aeq) || ~isempty(Beq) || ~isempty(Aineq)  || ~isempty(Bineq)
        subtype = 'linearconstraints';
    elseif ~isempty(LB) || ~isempty(UB)
        subtype = 'boundconstraints';
    else
        subtype = 'unconstrained';
    end
    % If Aeq or Aineq is not empty, then problem has linear constraints.
elseif ~isempty(Aeq) || ~isempty(Beq) || ~isempty(Aineq)  || ~isempty(Bineq)
    type = 'linearconstraints';
    % This condition satisfies bound constraints
elseif ~isempty(LB) || ~isempty(UB)
    type = 'boundconstraints';
    % If all constraints fields are empty then it is unconstrained
else
    type = 'unconstrained';
end

% Initialize output structure
OUTPUT = struct('function',FUN,'problemtype',type,'pollmethod',[], ...
    'searchmethod',[],'iterations',0,'funccount',0,'meshsize',[]);

% Store the random state
dflt = RandStream.getDefaultStream;
output.rngstate = struct('state',{dflt.State}, 'type',{dflt.Type});

% If nonlinear constraints, then subtype is needed to process linear
% constraints (see function preProcessLinearConstr)
if strcmp(type,'nonlinearconstr')
    type = subtype;
end

% Use default options if empty
if ~isempty(options) && ~isa(options,'struct')
    error('gads:patternsearch:optionsNotAStruct','Tenth input argument must be a valid structure created with PSOPTIMSET.');
elseif isempty(options)
    options = defaultopt;
end
% Save all user options
user_options = options;
% Validate all options
options = checkoptions(options,defaultopt,numberOfVariables);

% Bound correction
[LB,UB,msg,EXITFLAG] = checkbound(LB,UB,numberOfVariables);
if EXITFLAG < 0
    FVAL     =    [];
    X(:)     =  Iterate.x;
    OUTPUT.message = msg;
    if options.Verbosity > 0
        fprintf('%s\n',msg)
    end
    return
end
% Linear constraints correction
[Iterate.x,Aineq,Bineq,Aeq,Beq,LB,UB,msg,EXITFLAG] = ...
    preProcessLinearConstr(Iterate.x,Aineq,Bineq,Aeq,Beq,LB,UB,numberOfVariables,type,options.Verbosity);
if EXITFLAG < 0
    FVAL     =    [];
    X(:)     =  Iterate.x;
    OUTPUT.message = msg;
    if options.Verbosity > 0
        fprintf('%s\n',msg)
    end
    return
end
% Number of constraints
nineqcstr = size(Aineq,1);
neqcstr   = size(Aeq,1);
ncstr     = nineqcstr + neqcstr;

% Check the objective and constraint functions; Evaluate at the start point
[Iterate,OUTPUT.funccount] = poptimfcnchk(FUN,nonlcon,initialX,Iterate, ...
    options.Vectorized,objFcnArg,conFcnArg);

% Print some diagnostic information if verbosity > 2
if options.Verbosity > 2
    psdiagnose(FUN,nonlcon,initialX,nineqcstr,neqcstr,ncstr,user_options);
end
% Call appropriate private solver
switch(OUTPUT.problemtype)
    case 'unconstrained'
        [X,FVAL,EXITFLAG,OUTPUT] = pfminunc(FUN,objFcnArg,initialX,numberOfVariables,Iterate, ...
            options,defaultopt,OUTPUT);
    case 'boundconstraints'
        [X,FVAL,EXITFLAG,OUTPUT] = pfminbnd(FUN,objFcnArg,initialX,numberOfVariables,Iterate, ...
            LB,UB,options,defaultopt,OUTPUT);
    case 'linearconstraints'
        [X,FVAL,EXITFLAG,OUTPUT] = pfminlcon(FUN,objFcnArg,initialX,numberOfVariables,Iterate, ...
            Aineq,Bineq,Aeq,Beq,LB,UB,options,defaultopt,OUTPUT);
    case 'nonlinearconstr'
        [X,FVAL,EXITFLAG,OUTPUT] = pfmincon(FUN,objFcnArg,initialX,numberOfVariables,Iterate, ...
            Aineq,Bineq,Aeq,Beq,LB,UB,nonlcon,conFcnArg,options,defaultopt,OUTPUT,subtype);
end