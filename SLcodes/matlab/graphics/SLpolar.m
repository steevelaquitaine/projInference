%SLpolar.m
%
%       $Id: SLpolar.m $
%     usage: draw polar
%        by: steeve laquitaine
%      date: 140528
%   purpose: plot polar
%
%varargin are
%
%   'circle'
%   'radialCircles'
%   'xSpokes'
%   'SpokesTickSteps'
%   'area'

function hpol = SLpolar(varargin)
%POLAR  Polar coordinate plot.
%   POLAR(THETA, RHO) makes a plot using polar coordinates of
%   the angle THETA, in radians, versus the radius RHO.
%   POLAR(THETA, RHO, S) uses the linestyle specified in string S.
%   See PLOT for a description of legal linestyles.
%
%   POLAR(AX, ...) plots into AX instead of GCA.
%
%   H = POLAR(...) returns a handle to the plotted object in H.
%
%   Example:
%      t = 0 : .01 : 2 * pi;
%      polar(t, sin(2 * t) .* cos(2 * t), '--r');
%
%   See also PLOT, LOGLOG, SEMILOGX, SEMILOGY.

%   Copyright 1984-2013 The MathWorks, Inc.

% Parse possible Axes input
[cax, args, nargs] = axescheck(varargin{:});

%get distance between spokes in steps of values of theta
if sum(strcmp(varargin,'SpokesTickSteps'))==1
    SpokesTickSteps=varargin{find(strcmp(varargin,'SpokesTickSteps'))+1};
end

if sum(strcmp(varargin,'facecolor'))==1
    facecolor=varargin{find(strcmp(varargin,'facecolor'))+1};
end

if sum(strcmp(varargin,'edgecolor'))==1
    edgecolor=varargin{find(strcmp(varargin,'edgecolor'))+1};
end

if sum(strcmp(varargin,'linestyle'))==1
    linestyle=varargin{find(strcmp(varargin,'linestyle'))+1};
else
    linestyle='none';
end



theta = args{1};
rho = args{2};

%make sure theta and rho are row vectors
if size(theta,1)>size(theta,2)
    theta=theta';
end
if size(rho,1)>size(rho,2)
    rho=rho';
end

if ischar(rho)
    line_style = rho;
    rho = theta;
    [mr, nr] = size(rho);
    if mr == 1
        theta = 1 : nr;
    else
        th = (1 : mr)';
        theta = th(:, ones(1, nr));
    end
else
    line_style = 'auto';
end
if ischar(theta) || ischar(rho)
    error(message('MATLAB:polar:InvalidInputType'));
end
if ~isequal(size(theta), size(rho))
    error(message('MATLAB:polar:InvalidInputDimensions'));
end

% get hold state
cax = newplot(cax);

next = lower(get(cax, 'NextPlot'));
hold_state = ishold(cax);

% get x-axis text color so grid is in same color
if graphicsversion(cax, 'handlegraphics')
    tc = get(cax, 'XColor');
else
    % get the axis gridColor
    axGridColor = get(cax,'GridColor');
    tc = axGridColor;
end
ls = get(cax, 'GridLineStyle');

% Hold on to current Text defaults, reset them to the
% Axes' font attributes so tick marks use them.
fAngle = get(cax, 'DefaultTextFontAngle');
fName = get(cax, 'DefaultTextFontName');
fSize = get(cax, 'DefaultTextFontSize');
fWeight = get(cax, 'DefaultTextFontWeight');
fUnits = get(cax, 'DefaultTextUnits');
set(cax, ...
    'DefaultTextFontAngle', get(cax, 'FontAngle'), ...
    'DefaultTextFontName', get(cax, 'FontName'), ...
    'DefaultTextFontSize', get(cax, 'FontSize'), ...
    'DefaultTextFontWeight', get(cax, 'FontWeight'), ...
    'DefaultTextUnits', 'data');

% only do grids if hold is off
if ~hold_state
    
    % make a radial grid
    hold(cax, 'on');
    % ensure that Inf values don't enter into the limit calculation.
    arho = abs(rho(:));
    maxrho = max(arho(arho ~= Inf));
    hhh = line([-maxrho, -maxrho, maxrho, maxrho], [-maxrho, maxrho, maxrho, -maxrho], 'Parent', cax);
    set(cax, 'DataAspectRatio', [1, 1, 1], 'PlotBoxAspectRatioMode', 'auto');
    %v=[get(cax, 'XLim') get(cax, 'YLim')];
    
    %scale the axis limits to the max of the data
    v=[-maxrho maxrho -maxrho maxrho];
    
    ticks = sum(get(cax, 'YTick') >= 0);
    delete(hhh);
    
    %check radial limits and ticks
    rmin = 0;
    rmax = v(4);
    rticks = max(ticks - 1, 2);
    if rticks > 5   % see if we can reduce the number
        if rem(rticks, 2) == 0
            rticks = rticks / 2;
        elseif rem(rticks, 3) == 0
            rticks = rticks / 3;
        end
    end
    
    %define a circle
    th = 0 : pi / 50 : 2 * pi;
    xunit = cos(th);
    yunit = sin(th);
    
    %now really force points on x/y axes to lie on them exactly
    inds = 1 : (length(th) - 1) / 4 : length(th);
    xunit(inds(2 : 2 : 4)) = zeros(2, 1);
    yunit(inds(1 : 2 : 5)) = zeros(3, 1);
    
    %plot background if necessary
    if strcmp(varargin,'circle')
        if ~ischar(get(cax, 'Color'))
            patch('XData', xunit * rmax, 'YData', yunit * rmax, ...
                'EdgeColor', tc, 'FaceColor', 'none', ...
                'HandleVisibility', 'off', 'Parent', cax);
        end
    end
    
    %(case we want to draw radial circles)
    if strcmp(varargin,'radialCircles')
        c82 = cos(82 * pi / 180);
        s82 = sin(82 * pi / 180);
        rinc = (rmax - rmin) / rticks;
        for i = (rmin + rinc) : rinc : rmax
            hhh = line(xunit * i, yunit * i, 'LineStyle', ls, 'Color', tc, 'LineWidth', 1, ...
                'HandleVisibility', 'on', 'Parent', cax);
            text((i + rinc / 20) * c82, (i + rinc / 20) * s82, ...
                ['  ' num2str(i)], 'VerticalAlignment', 'bottom', ...
                'HandleVisibility', 'on', 'Parent', cax);
        end
        
        %Make outer circle solid
        set(hhh, 'LineStyle', '-');
    end
    
    %plot spokes
    if sum(strcmp(varargin,'xSpokes'))==1
        
        %get a bit of space to data
        rmax=1.1*rmax;
        
        %get circular data in radian
        th=theta;
        
        %convert to cartesian coordinates (something wrong here)
        %x coordinate
        cst=cos(th);
        
        %y coordinate
        snt=sin(th);

        %set directional vector's origin coordinates [0,0]
        xo=zeros(1,numel(th));
        yo=zeros(1,numel(th));
        
        %set directional vector's end coordinates [0,0]
        line([xo;rmax*cst], [yo;rmax*snt], 'LineStyle', ls, 'Color', [.7 .7 .7], 'LineWidth', 1, ...
            'HandleVisibility', 'on', 'Parent', cax);
        
        %annotate spokes in degrees
        rt=1.1*rmax;
        for i=1:SpokesTickSteps:length(th)
            text(rt*cst(i),rt*snt(i),int2str(SLra2d(th(i))),...
                'HorizontalAlignment', 'center', ...
                'HandleVisibility', 'on', 'Parent', cax,...
                'fontsize',5);
            if i == length(th)
                loc = int2str(0);
            else
                loc = int2str(180 + SLra2d(th));
            end
        end
        
        %default spokes 
    else
        %th=(1 : 6) * 2 * pi / 12;
        th=deg2rad([0 90 180 270]);
        cst = cos(th);
        snt = sin(th);
        cs = [-cst; cst];
        sn = [-snt; snt];
        line(rmax * cs, rmax * sn, 'LineStyle', ls, 'Color', [.7 .7 .7], 'LineWidth', 1, ...
            'HandleVisibility', 'on', 'Parent', cax);
       
        %annotate spokes in degrees
        rt=1.1*rmax;
        for i=1:length(th)
            text(rt*cst(i),rt*snt(i),int2str(SLra2d(th(i))),...
                'HorizontalAlignment', 'center', ...
                'HandleVisibility', 'on', 'Parent', cax);
            %if i == length(th)
            %    loc = int2str(0);
            %else
            %    loc = int2str(180 + SLra2d(th));
            %end
            %text(-rt * cst(i), -rt * snt(i), loc, 'HorizontalAlignment', 'center', ...
            %    'HandleVisibility', 'on', 'Parent', cax);
        end
    end
    
    %set view to 2-D
    view(cax,2);
    
    % set axis limits
    axis(cax,rmax*[-1, 1, -1.15, 1.15]);
end

% Reset defaults.
set(cax, ...
    'DefaultTextFontAngle', fAngle , ...
    'DefaultTextFontName', fName , ...
    'DefaultTextFontSize', fSize, ...
    'DefaultTextFontWeight', fWeight, ...
    'DefaultTextUnits', fUnits );

% transform data to Cartesian coordinates.
xx = rho .* cos(theta);
yy = rho .* sin(theta);

%case area plot
if sum(strcmp(varargin,'area'))==1
    q=area(xx, yy, 'Parent', cax,'FaceColor',facecolor,...
        'lineStyle','none');
    
    q2=plot(xx, yy, 'Parent', cax,'color',facecolor,...
        'linesmoothing','on','linewidth',1,...
        'lineStyle',linestyle);
    
    %plot data on top of grid as default
else
    if strcmp(line_style, 'auto')
        q = plot(xx, yy, 'Parent', cax,'color',edgecolor,...
            'linesmoothing','on','linewidth',2,'linestyle',linestyle);
    else
        q = plot(xx, yy, line_style, 'Parent', cax,...
            'color','k','linesmoothing','on','linewidth',2,'lineStyle',linestyle);
    end
end

if nargout == 1
    hpol = q;
end

if ~hold_state
    set(cax, 'DataAspectRatio', [1, 1, 1]), axis(cax, 'off');
    set(cax, 'NextPlot', next);
end
set(get(cax, 'XLabel'), 'Visible', 'on');
set(get(cax, 'YLabel'), 'Visible', 'on');

if ~isempty(q) && ~isdeployed
    makemcode('RegisterHandle', cax, 'IgnoreHandle', q,...
        'FunctionName', 'polar');
end


