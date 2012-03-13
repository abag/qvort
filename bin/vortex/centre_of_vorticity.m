%read in the cov_info (centre of vorticity) file and plot various cov information
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
A=load('./data/centre_of_vorticity_info.log');
t=A(:,1) ; covx=A(:,2) ; covy=A(:,3) ; covz=A(:,4) ;
covux=A(:,5) ; covuy=A(:,6) ; covuz=A(:,7) ;
switch option
  case 'print'
    figure('visible','off');
  otherwise
    figure('Name', 'Centre of vorticity - position')      
end
  subplot(3,1,1)
    plot(t,covx,'-r','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('x','FontSize',14)
  subplot(3,1,2)
    plot(t,covy,'-b','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('y','FontSize',14)
  subplot(3,1,3)
    plot(t,covz,'-m','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('z','FontSize',14)
if option=='print'
    disp('printing to cov_position_info.eps')
    print('-depsc','./cov_position_info.eps')
end
switch option
  case 'print'
    figure('visible','off');
  otherwise
    figure('Name', 'Centre of vorticity - velocity')      
end
  subplot(3,1,1)
    plot(t,covux,'-r','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('u(x)','FontSize',14)
  subplot(3,1,2)
    plot(t,covuy,'-b','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('u(y)','FontSize',14)
  subplot(3,1,3)
    plot(t,covuz,'-m','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('u(z)','FontSize',14)
if option=='print'
    disp('printing to cov_velocity_info.eps')
    print('-depsc','./cov_velocity_info.eps')
end

