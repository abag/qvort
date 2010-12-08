function particle_adpative(option)
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
A=load('data/qp_eta.log');
t=A(:,1) ; eta=A(:,2) ; dt=A(:,3) ;
switch option
  case 'print'
    figure('visible','off');
  otherwise
    figure('Name', 'adaptive details')
end
  subplot(2,1,1)
    semilogy(t,dt,'-','LineWidth',1,'Color',rgb('Indigo'));
    set(gca,'FontSize',14);
    xlabel('t','FontSize',14);
    ylabel('log(dt)','FontSize',14);
  subplot(2,1,2)
    plot(t,eta,'-','LineWidth',1,'Color','r');
    set(gca,'FontSize',14);
    xlabel('t','FontSize',14);
    ylabel('max(\theta_{quasi})','FontSize',14);
switch option
  case 'print'
    disp('printing to qp_adaptive.eps')
    print('-depsc','./qp_adaptive.eps')
end
