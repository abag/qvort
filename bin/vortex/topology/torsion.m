%script to read in the torsion information and plot
%if given the option print will print to .eps file rather than screen
function torsion(option)
if nargin==0     
  option='empty';
end
switch option
case 'print'
    disp('will not print to screen but instead to .eps files')
case 'empty'
    otherwise
    disp('incorrect option, aborting script and printing help:')
    help torsion
    return
end
A=load('./data/torsion.log');
B=load('./data/ts.log') ; 
t=B(:,2) ; tbar=A(:,1) ; tmin=A(:,2) ; tmax=A(:,3) ;
switch option
  case 'print'
    figure('visible','off');
  otherwise
    figure('Name', 'torsion')      
end
  subplot(3,1,1)
    plot(t,tbar,'-r','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('mean tors','FontSize',14)
  subplot(3,1,2)
    plot(t,tmin,'-b','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('min tors','FontSize',14)
  subplot(3,1,3)
    plot(t,tmax,'-m','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('max tors','FontSize',14)
switch option
  case 'print'
    disp('printing to torsion.eps')
    print('-depsc','./torsion.eps')
end
