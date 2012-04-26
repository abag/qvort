%read in the ts file and plot reconnection information
%you must have set recon_info T in run.in
%if given the option print will print to .eps file rather than screen
function recon_rate(option)
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
A=load('./data/ts.log');
t=A(:,2) ; rcount=A(:,4) ; rmcount=A(:,11) ;
recon_rate=gradient(rcount,t(2)-t(1));
rm_rate=gradient(rmcount,t(2)-t(1));
switch option
  case 'print'
    figure('visible','off');
  otherwise
    figure('Name', 'reconnection rates')      
end
  subplot(2,2,1)
    plot(t,recon_rate,'-r','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('recon rate','FontSize',14)
  subplot(2,2,2)
    plot(t,rcount,'-r','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('recon count','FontSize',14)
  subplot(2,2,3)
    plot(t,rm_rate,'-b','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('rm rate','FontSize',14)
  subplot(2,2,4)
    plot(t,rmcount,'-b','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('rm count','FontSize',14)
if option=='print'
    disp('printing to recon_rate.eps')
    print('-depsc','./recon_rate.eps')
end
switch option
  case 'print'
    figure('visible','off');
  otherwise
    figure('Name', 'recon ratios')      
end
  subplot(2,1,1)
    plot(t,recon_rate./rm_rate,'-b','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('recon/rm rate','FontSize',14)
  subplot(2,1,2)
    plot(t,rcount./rmcount,'-r','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('recon/rm count','FontSize',14)
if option=='print'
  disp('printing to recon_ratios.eps')
  print('-depsc','./recon_ratios.eps')
end
