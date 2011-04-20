%read in the particle ts file and plot various dignostic information
%if given the option print will print to .eps file rather than screen
function par_ts(option)
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
B=load('data/par_ts.log');
t=B(:,2) ; pmaxu=B(:,3) ; pmaxdu=B(:,4) ; purms=B(:,5) ; psep=B(:,6) ;
switch option
  case 'print'
    figure('visible','off');
  otherwise
    figure('Name', 'particle information')      
end
subplot(2,2,1)
 plot(t,pmaxu,'-c','LineWidth',2);
 set(gca,'FontSize',14)
   xlabel('t','FontSize',14)
   ylabel('max u','FontSize',14)
subplot(2,2,2)
 plot(t,pmaxdu,'-r','LineWidth',2);
   set(gca,'FontSize',14)
   xlabel('t','FontSize',14)
   ylabel('max du/dt','FontSize',14)
subplot(2,2,3)
 plot(t,purms,'-y','LineWidth',2);
   set(gca,'FontSize',14)
   xlabel('t','FontSize',14)
   ylabel('urms','FontSize',14)
subplot(2,2,4)
 plot(t,psep,'-g','LineWidth',2);
   set(gca,'FontSize',14)
   xlabel('t','FontSize',14)
   ylabel('particle sep.','FontSize',14)
if option=='print'
  disp('printing to particle_information.eps')
  print('-depsc','./particle_information.eps')
end
 
