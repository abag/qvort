%read in the B_ts file and plot various dignostic information
%if given the option print will print to .eps file rather than screen
function B_ts(option)
if nargin==0     
  option='empty';
end
switch option
case 'print'
  disp('will not print to screen but instead to .eps files')
case 'empty'
  otherwise
  disp('incorrect option, aborting script and printing help:')
  help ts
  return
end
A=load('data/B_ts.log');
t=A(:,1) ; B_max=A(:,2) ; B_min=A(:,3) ; B_rms=A(:,4) ; B_E=A(:,5); 
switch option
  case 'print'
    figure('visible','off');
  otherwise
    figure('Name', 'magnetic information')      
end
  subplot(2,2,1)
    plot(t,B_max,'-r','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('B_{max}','FontSize',14)
  subplot(2,2,2)
    plot(t,B_min,'-b','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('B_{min}','FontSize',14)
  subplot(2,2,3)
    plot(t,B_rms,'-m','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('B_{rms}','FontSize',14)
  subplot(2,2,4)
    plot(t,B_E,'-k','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('E_{mag}','FontSize',14)
if option=='print'
    disp('printing to magnetic_information.eps')
    print('-depsc','./magnetic_information.eps')
end
