%read in the ts file and plot
A=load('data/ts.log');
t=A(:,2) ; pcount=A(:,3) ; rcount=A(:,4) ; sep=A(:,5) ; l=A(:,6) ; 
maxu=A(:,7) ; maxdu=A(:,8) ;
figure('Name', 'filament information')
  subplot(2,2,1)
    plot(t,pcount,'-r','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('pcount','FontSize',14)
  subplot(2,2,2)
    plot(t,rcount,'-b','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('recon count','FontSize',14)
  subplot(2,2,3)
    plot(t,sep,'-m','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('avg sep','FontSize',14)
  subplot(2,2,4)
    plot(t,l,'-g','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('L','FontSize',14)
figure('Name', 'velocity information')
  subplot(2,1,1)
    plot(t,maxu,'-b','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('max(u)','FontSize',14)
  subplot(2,1,2)
    plot(t,maxdu,'-r','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('max(du)','FontSize',14)
 if exist('data/par_ts.log');
   B=load('data/par_ts.log');
   t=B(:,2) ; pmaxu=B(:,3) ; pmaxdu=B(:,4) ; purms=B(:,5) ; psep=B(:,6) ;
   figure('Name', 'particle information')
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
 end
