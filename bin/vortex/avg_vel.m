%read in the basic_velocity_info file and plot various velocity information
%if given the option print will print to .eps file rather than screen
function avg_vel(option)
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
A=load('./data/basic_velocity_info.log');
t=A(:,1) ; avgux=A(:,5) ; avguy=A(:,6) ; avguz=A(:,7) ; avgu=A(:,3) ;

switch option
  case 'print'
    figure('visible','off');
  otherwise
    figure('Name', 'Basic velocity information')      
end
  subplot(4,1,1)
    plot(t,avgux,'-r','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('avg(u_x)','FontSize',14)
  subplot(4,1,2)
    plot(t,avguy,'-b','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('avg(u_y)','FontSize',14)
  subplot(4,1,3)
    plot(t,avguz,'-m','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('avg(u_z)','FontSize',14)
  subplot(4,1,4)
    plot(t,avgu,'-m','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('avg(u)','FontSize',14)
if option=='print'
    disp('printing to basic_velocity_info.eps')
    print('-depsc','./basic_velocity_info.eps')
end
