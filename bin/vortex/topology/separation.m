%script to read in the separation information and plot
%if given the option print will print to .eps file rather than screen
function separation(option)
dims=load('./data/dims.log');
delta=dims(1);
if nargin==0     
  option='empty';
end
switch option
case 'print'
    disp('will not print to screen but instead to .eps files')
case 'empty'
    otherwise
    disp('incorrect option, aborting script and printing help:')
    help separation
    return
end
A=load('data/sep_info.log');
t=A(:,1) ; sbar=A(:,2) ; smin=A(:,3) ; smax=A(:,4) ;
switch option
  case 'print'
    figure('visible','off');
  otherwise
    figure('Name', 'separation')      
end
  subplot(3,1,1)
    plot(t,sbar,'-r','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('mean sep','FontSize',14)
  subplot(3,1,2)
    plot(t,smin,'-b','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('min sep','FontSize',14)
  subplot(3,1,3)
    plot(t,smax,'-k','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('max sep','FontSize',14)
switch option
  case 'print'
    disp('printing to separation.eps')
    print('-depsc','./separation.eps')
end
