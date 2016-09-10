function [h_main, h_inset]=inset(main_handle, inset_handle,inset_size)

%The function plotting figure inside figure (main and inset) from 2 existing figures.
%inset_size is the fraction of inset-figure size, default value is 0.35
%The outputs are the axes-handles of both.
% 
%An examle can found in the file: inset_example.m
%close all
%x1=0:50;
%y1=besselj(0,x1);
%x2=2.35:.001:2.45;
%y2=besselj(0,x2);
%fig1=figure(1);
%plot(x1,y1,'b','linewidth',2)
%hold on
%plot(x1,0*x1,':k')
%set(gca,'fontsize',15)
%title ('bessel function')
%xlabel('X')
%ylabel('Y')
%fig2=figure(2);
%plot(x2,y2,'r')
%hold on
%plot(x2,0*x2,':k')
%title ('close-up')
%[h_m h_i]=inset(fig1,fig2);

%set(h_i,'xtick',2.35:.025:2.45,'xlim',[2.35,2.45])
% The script is an example for the function "inset"
% figure 1 is a plot of 'zero order bessel function of the first kind' in the range of x=0 to x=50.
% figure 2 is a close-up of this function around its first zero.
% figure 3 shows the function in the main plot, and the close-up in the
% inset plot.
%
% Moshe Lindner, August 2010 (C).



if nargin==2
    inset_size=0.35;
end

inset_size=inset_size*.7;
figure
new_fig=gcf;
main_fig = findobj(main_handle,'Type','axes');
h_main = copyobj(main_fig,new_fig);
set(h_main,'Position',get(main_fig,'Position'))
inset_fig = findobj(inset_handle,'Type','axes');
h_inset = copyobj(inset_fig,new_fig);
ax=get(main_fig,'Position');
set(h_inset,'Position', [.7*ax(1)+ax(3)-inset_size .5*ax(2)+ax(4)-inset_size inset_size inset_size])