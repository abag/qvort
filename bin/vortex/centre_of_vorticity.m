%read in the com_info (centre of mass) file and plot various com information
%if given the option print will print to .eps file rather than screen
function centre_of_vorticity(option)
if nargin==0     
  option='empty';
end
switch option
case 'print'
  disp('will not print to screen but instead to .eps files')
case 'empty'
  otherwise
  disp('incorrect option, aborting script')
  return
end
A=load('./data/daniel_com_info.log');
t=A(:,1) ; comx=A(:,2) ; comy=A(:,3) ; comz=A(:,4) ;

switch option
  case 'print'
    figure('visible','off');
  otherwise
    figure('Name', 'Centre of mass information')      
end
  subplot(3,1,1)
    plot(t,comx,'-r','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('com(x)','FontSize',14)
  subplot(3,1,2)
    plot(t,comy,'-b','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('com(y)','FontSize',14)
  subplot(3,1,3)
    plot(t,comz,'-m','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('com(z)','FontSize',14)
if option=='print'
    disp('printing to com_info.eps')
    print('-depsc','./com_info.eps')
end
