function particle_energy(option)
A=load('data/qp_energy.log');
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
switch option
  case 'print'
    figure('visible','off');
  otherwise
    figure('Name', 'quasi-particle energy')
end
plot(A(:,1),A(:,2),'LineWidth',2)
  set(gca,'FontSize',14)
  xlabel('t','FontSize',14)
  ylabel('energy','FontSize',14)
switch option
  case 'print'
    disp('printing to qp_energy.eps')
    print('-depsc','./qp_energy.eps')
end




