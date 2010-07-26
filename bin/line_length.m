%read in the ts file and plot
A=load('data/ts.log');
t=A(:,2) ; l=A(:,6) ; 
scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(3)/1.2 scrsz(4)/2],'PaperPosition',[0.25 2.5 28.0 12.0],'color','w','visible','on','Name', 'line_length')
  subplot(1,2,1)
    plot(t,l,'-r','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('l','FontSize',14)
  subplot(1,2,2)
    p=polyfit(t,log(l),1);
    dummy_ll=p(1)*t+p(2);
    plot(t,log(l),'-b',t,dummy_ll,'-.r','LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('t','FontSize',14)
    ylabel('log(l)','FontSize',14)
