function [ output_args ] = MakePretty( input_args )
%MAKEPUBLISHABLE Summary of this function goes here
%   Detailed explanation goes here


% Set figure's dimensions
% --------------------------------------------------------------------
% we set the units of the measures used through the file
%
% [ inches | centimeters | normalized | points | {pixels} | characters ]
set(gcf, 'Units', 'centimeters');
% we set the position and dimension of the figure ON THE SCREEN
%
% NOTE: measurement units refer to the previous settings!
afFigurePosition = [0 0 25 25];         % [pos_x pos_y width_x width_y]
set(gcf, 'Position', afFigurePosition);  % [left bottom width height]
% we link the dimension of the figure ON THE PAPER in such a way that
% it is equal to the dimension on the screen
%
% ATTENTION: if PaperPositionMode is not 'auto' the saved file
% could have different dimensions from the one shown on the screen!
set(gcf, 'PaperPositionMode', 'auto');
% in order to make matlab to do not "cut" latex-interpreted axes labels
set(gca, 'Units','normalized',... %
 'Position',[0.15 0.2 0.75 0.7]);


% Get rid of white margin
set(gca,'LooseInset',get(gca,'TightInset'))


% here we select which output file extension we want
bPrintOnFile_Pdf    = 1; % [0 (false)
bPrintOnFile_Eps    = 0; % [0 (false)
% we select the file path
%
% NOTE: do NOT insert extensions!
%strFilePath = '../images/my_figure';
strFilePath = 'figures';
% we select the printing resolution
iResolution = 600;

% we select to crop or not the figure
bCropTheFigure = 1; % [0 (false)   1 (true)]
% ATTENTION: if PaperPositionMode is not 'auto' the saved file
% could have different dimensions from the one shown on the screen!
set(gcf, 'PaperPositionMode', 'auto');
% saving on file: requires some checks
if( bPrintOnFile_Pdf || bPrintOnFile_Eps )
    %
    % NOTE: if you want a .pdf with encapsulated fonts you need to save an
    % .eps and then convert it => it is always necessary to produce the .eps
    %
    % if we want to crop the figure we do it
    if( bCropTheFigure )
        print('-depsc2', sprintf('-r%d', iResolution), strcat(strFilePath, '.eps'));
    else
        print('-depsc2', '-loose', sprintf('-r%d', iResolution), strcat(strFilePath, '.eps'));
    end;
%

% if we want the .pdf we produce it
if( bPrintOnFile_Pdf )
    %
    % here we convert the .eps encapsulating the fonts
    system(...
          sprintf(...
            'epstopdf --gsopt=-dPDFSETTINGS=/prepress --outfile=%s.pdf  %s.eps',   ...
            strFilePath,                                                           ...
            strFilePath));
% end;
    %
    % if we do not want the .eps we remove it
    if( ~bPrintOnFile_Eps )
        delete(sprintf('%s.eps', strFilePath));
    end;
    %
end;% saving on file



end

