%script to read in the curvature information and plot
%if given the option print will print to .eps file rather than screen
function curvature(option)
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
    help curvature
    return
end
A=load('./data/curvature.log');
B=load('./data/ts.log') ; 
t=B(:,2) ; cbar=A(:,1) ; cmin=A(:,2) ; cmax=A(:,3) ;
cmaxmax(1:length(A))=sqrt(3)/delta;
%cmaxmax(1:length(A))=max(cmax)-20
switch option
  case 'print'
    figure('visible','off');
  otherwise
    figure('Name', 'curvature')      
end
  subplot(3,1,1)
    plot(t,cbar,'-r','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('mean curv','FontSize',14)
  subplot(3,1,2)
    plot(t,cmin,'-b','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('min curv','FontSize',14)
  subplot(3,1,3)
    plot(t,cmax,'-m',t,cmaxmax,'--k','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('max curv','FontSize',14)
switch option
  case 'print'
    disp('printing to curvature.eps')
    print('-depsc','./curvature.eps')
end
