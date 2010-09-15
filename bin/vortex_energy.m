%script to read in the vortex energy information and plot
%if given the option print will print to .eps file rather than screen
function vortex_energy(option)
if nargin==0     
  option='empty';
end
switch option
case 'print'
    disp('will not print to screen but instead to .eps files')
case 'empty'
    otherwise
    disp('incorrect option, aborting script and printing help:')
    help vortex_energy
    return
end
A=load('data/energy.log');
B=load('data/ts.log') ; 
t=B(:,2) ; energy=A(:,1) ; 
switch option
  case 'print'
    figure('visible','off');
  otherwise
    figure('Name', 'vortex energy')      
end
  plot(t,energy,'-r','LineWidth',2);
  set(gca,'FontSize',14)
  xlabel('t','FontSize',14)
  ylabel('energy','FontSize',14)
switch option
  case 'print'
    disp('printing to vortex_energy.eps')
    print('-depsc','./vortex_energy.eps')
end