%script to read in the line separation information and plot
%if given the option print will print to .eps file rather than screen
function line_separation(option)
if nargin==0     
  option='empty';
end
switch option
case 'print'
    disp('will not print to screen but instead to .eps files')
case 'empty'
    otherwise
    disp('incorrect option, aborting script and printing help:')
    help curvature
    return
end
A=load('./data/line_sep_info.log');
t=A(:,1) ; sep_bar=A(:,2) ; sep_max=A(:,3) ; sep_min=A(:,4) ;
switch option
  case 'print'
    figure('visible','off');
  otherwise
    figure('Name', 'line separation')      
end
  subplot(3,1,1)
    plot(t,sep_bar,'-r','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('mean sep','FontSize',14)
  subplot(3,1,2)
    plot(t,sep_min,'-b','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('min sep','FontSize',14)
  subplot(3,1,3)
    plot(t,sep_max,'-m','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('max sep','FontSize',14)
switch option
  case 'print'
    disp('printing to line_sep.eps')
    print('-depsc','./line_sep.eps')
end
