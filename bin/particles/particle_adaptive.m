A=load('data/qp_eta.log');
t=A(:,1) ; eta=A(:,2) ; dt=A(:,3) ;
figure('Name', 'adaptive details')
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