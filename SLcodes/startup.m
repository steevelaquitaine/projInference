
%steeve laquitaine

%add my codes path
% addpath(genpath('~/Dropbox/myDropbox/Codes'))
% addpath(genpath('~/proj/'))

% %matlab theme
% com.mathworks.services.Prefs.setBooleanPref('ColorsUseSystem',0);
% % com.mathworks.services.Prefs.setColorPref('ColorsBackground',java.awt.Color(0,0,0));
% com.mathworks.services.Prefs.setColorPref('ColorsText',java.awt.Color(133/255,153/255,0/255));
% com.mathworks.services.Prefs.setColorPref('Colors_M_Comments',java.awt.Color(181/255,137/255,0/255));
% com.mathworks.services.Prefs.setColorPref('Colors_M_Keywords',java.awt.Color(38/255,139/255,210/255));
% com.mathworks.services.Prefs.setColorPref('Colors_M_Strings',java.awt.Color(42/255,161/255,152/255));
% 
% % % Set Matlab color preferences to the 'classic' Emacs color scheme
% % % Setup Color Scheme
% %     sycb = 0;               %Set 'Use system colors' checkbox to False
% txc  = [  1   1  1  ];  %Set 'Text' color
% bgc  = [0.1 0.1 0.1];   %Set 'Background' color
% %     kwd  = [  0  1   1  ];  %Set 'Keywords' color
% cmt  = [ .5 .5  .5 ];  %Set 'Comments' color
% %     str  = [ 0.804 0.506 0.384 ];  % Set 'Strings' color
% %     ustr = [ 0.804 0.506 0.384 ];  % Set 'Unterminated strings' color
% %     scmd = [ 0 0 1.000 ];  % Set 'System commands' color
% %     errs = [ 0.800 0.310 0.263 ];  % Set 'Errors' color
% %     hyp  = [ 0.298 0.678 1.000 ];  % Set 'Hyperlinks' color
% %     warn = [ 1.000 0.627 0.478 ];  % Set 'Warnings' color
% %     afhb = .7;           % Set 'Autofix highlight' checkbox to True
% %     afh  = [ 0.584 0.745 0.847 ];  % Set 'Autofix highlight' color
% %     ahib = 1;           % Set 'Automatically highlight' checkbox to False
% %     ahi  = [ 0.506 0.588 0.604 ];  % Set 'Automatically highlight' color
% %     vwss = [ 1.000 0.627 0.478 ];   % Set 'Variables with shared scope' color
% %     clhb = 0;           % Set 'Highlight cells' checkbox to False
% %     hclb = 0;           % Set 'Highlight current line' checkbox to False
% %     hcl  = [ 0.506 0.588 0.604 ];  % Set 'Highlight current line' color
% %     slnb = 1;           % Set 'Show line numbers' checkbox to True
% %     rtlb = 0;           % Set 'Show line' checkbox in Right-hand text limit to False
% %
% % %Reset colors
% % %     sycb = 1;           % Set 'Use system colors' checkbox to True
% % % Need to get colors setup to return to default
% %
% % %Desktop tool colors
% % com.mathworks.services.Prefs.setBooleanPref('ColorsUseSystem',sycb); clear('sycb')
% %
% com.mathworks.services.Prefs.setColorPref('ColorsText',java.awt.Color(txc(1), txc(2), txc(3)));
% com.mathworks.services.ColorPrefs.notifyColorListeners('ColorsText'); clear('txc')
% 
% com.mathworks.services.Prefs.setColorPref('ColorsBackground',java.awt.Color(bgc(1), bgc(2), bgc(3)));
% com.mathworks.services.ColorPrefs.notifyColorListeners('ColorsBackground'); clear('bgc')
% %
% % %Setup MATLAB syntax highlighting colors
% % com.mathworks.services.Prefs.setColorPref('Colors_M_Keywords',java.awt.Color(kwd(1), kwd(2), kwd(3)));
% % com.mathworks.services.ColorPrefs.notifyColorListeners('Colors_M_Keywords'); clear('kwd')
% %
% com.mathworks.services.Prefs.setColorPref('Colors_M_Comments',java.awt.Color(cmt(1), cmt(2), cmt(3)));
% com.mathworks.services.ColorPrefs.notifyColorListeners('Colors_M_Comments'); clear('cmt')
% %
% % com.mathworks.services.Prefs.setColorPref('Colors_M_Strings',java.awt.Color(str(1), str(2), str(3)));
% % com.mathworks.services.ColorPrefs.notifyColorListeners('Colors_M_Strings'); clear('str')
% %
% % com.mathworks.services.Prefs.setColorPref('Colors_M_UnterminatedStrings',java.awt.Color(ustr(1), ustr(2), ustr(3)));
% % com.mathworks.services.ColorPrefs.notifyColorListeners('Colors_M_UnterminatedStrings'); clear('ustr')
% %
% % com.mathworks.services.Prefs.setColorPref('Colors_M_SystemCommands',java.awt.Color(scmd(1), scmd(2), scmd(3)));
% % com.mathworks.services.ColorPrefs.notifyColorListeners('Colors_M_SystemCommands'); clear('scmd')
% %
% % com.mathworks.services.Prefs.setColorPref('Colors_M_Errors',java.awt.Color(errs(1), errs(2), errs(3)));
% % com.mathworks.services.ColorPrefs.notifyColorListeners('Colors_M_Errors');clear('errs')
% %
% % % Other colors
% % com.mathworks.services.Prefs.setColorPref('Colors_HTML_HTMLLinks',java.awt.Color(hyp(1), hyp(2), hyp(3)));
% % com.mathworks.services.ColorPrefs.notifyColorListeners('Colors_HTML_HTMLLinks'); clear('hyp')
% %
% % % Code analyzer colors
% % com.mathworks.services.Prefs.setColorPref('Colors_M_Warnings',java.awt.Color(warn(1), warn(2), warn(3)));
% % com.mathworks.services.ColorPrefs.notifyColorListeners('Colors_M_Warnings'); clear('warn')
% %
% % com.mathworks.services.Prefs.setBooleanPref('ColorsUseMLintAutoFixBackground',afhb);
% % com.mathworks.services.Prefs.setColorPref('ColorsMLintAutoFixBackground',java.awt.Color(afh(1), afh(2), afh(3)));
% % com.mathworks.services.ColorPrefs.notifyColorListeners('ColorsMLintAutoFixBackground'); clear('afhb','afh')
% %
% % % Variable and function colors
% % com.mathworks.services.Prefs.setBooleanPref('Editor.VariableHighlighting.Automatic',ahib);
% % com.mathworks.services.Prefs.setColorPref('Editor.VariableHighlighting.Color',java.awt.Color(ahi(1), ahi(2), ahi(3)));
% % com.mathworks.services.ColorPrefs.notifyColorListeners('Editor.VariableHighlighting.Color'); clear('ahib','ahi')
% %
% % % *** need to add checkbox option for variables with shared scope here
% % com.mathworks.services.Prefs.setColorPref('Editor.NonlocalVariableHighlighting.TextColor',java.awt.Color(vwss(1), vwss(2), vwss(3)));
% % com.mathworks.services.ColorPrefs.notifyColorListeners('Editor.NonlocalVariableHighlighting.TextColor');clear('vwss')
% %
% % % Cell display options
% % com.mathworks.services.Prefs.setBooleanPref('EditorCodepadHighVisible',clhb); clear('clhb')
% % % Highlight cells color is 'Editorhighlight-lines'
% %
% % % Editor/Debugger General display options
% % com.mathworks.services.Prefs.setBooleanPref('Editorhighlight-caret-row-boolean',hclb);
% % com.mathworks.services.Prefs.setColorPref('Editorhighlight-caret-row-boolean-color',java.awt.Color(hcl(1), hcl(2), hcl(3)));
% % com.mathworks.services.ColorPrefs.notifyColorListeners('Editorhighlight-caret-row-boolean-color'); clear('hclb','hcl')
% % com.mathworks.services.Prefs.setBooleanPref('EditorShowLineNumbers',slnb); clear('slnb')
% % com.mathworks.services.Prefs.setBooleanPref('EditorRightTextLineVisible',rtlb);
% % %EditorRightTextLineLimit=I80
% % %EditorRightTextLimitLineWidth=I1
% 
% % Cleanup
% clear('sol'); clear('rtlb')