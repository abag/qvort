%script to read in the loop counter information and plot
%if given the option print will print to .eps file rather than screen
function loop_counter(option)
if nargin==0     
  option='empty';
end
switch option
case 'print'
    disp('will not print to screen but instead to .eps files')
case 'empty'
    otherwise
    disp('incorrect option, aborting script and printing help:')
    help loop_counter
    return
end
A=load('./data/loop_counter.log');
t=A(:,1) ; 
nloops=A(:,2) ; 
switch option
  case 'print'
    figure('visible','off');
  otherwise
    figure('Name', 'loop counter')      
end
  plot(t,nloops,'-b','LineWidth',2);
  set(gca,'FontSize',14)
  xlabel('t','FontSize',14)
  ylabel('N (loops)','FontSize',14)
switch option
  case 'print'
    disp('printing to loop_counter.eps')
    print('-depsc','./loop_counter.eps')
end
