function CF_compare(runcount)
for i=1:runcount
    fIN=sprintf('./CF_Tree_spec%1d/data/ts.log',i);
    b=load(fIN);
    t=b(:,2);
    L=b(:,6);
    nrecon=b(:,4);
    recon_rate=gradient(nrecon,t(2)-t(1));
    if i==1
      h1=figure('Name','line length');
    else
      figure(h1)
      hold on
    end
    plot(t,L)
    if i==1
      h2=figure('Name','recon rate');
    else
      figure(h2)
      hold on
    end
    semilogy(t,recon_rate)
    %plot(t,recon_rate)
    clear b
end
figure(h1)
xlabel('t','FontSize',16)
ylabel('line length','FontSize',16)
set(gca,'FontSize',16)
hold off
figure(h2)
xlabel('t','FontSize',16)
ylabel('recon rate','FontSize',16)
set(gca,'FontSize',16)
hold off
