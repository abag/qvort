function CF_compare(runcount)
colormap(hsv(6));
counter=1;
for i=1:runcount
    if i==4 || i==5
        continue
    end
    fIN=sprintf('./CF_Tree_spec%1d/data/ts.log',i);
    fIN_L=sprintf('./CF_Tree_spec%1d/LxLyLz',i);
    fIN_A=sprintf('./CF_Tree_spec%1d/Anisotropy',i);
    b=load(fIN);
    load(fIN_L);
    load(fIN_A);
    t=b(:,2);
    L=b(:,6);
    mean_curv=b(:,10);
    Density=L/0.001;
    nrecon=b(:,4);
    recon_rate=gradient(nrecon,t(2)-t(1));
    if i==1
      h1=figure('Name','line density');
    else
      figure(h1)
      hold on
    end
    %plot(t,Density,'LineWidth',1.5,'Color',rgb('HotPink'));
    if i==1
        plot(t,Density,'LineWidth',1.5,'Color',rgb('LimeGreen'));
    elseif i==2
        plot(t,Density,'LineWidth',1.5,'Color',rgb('Orange'));
    elseif i==3
        plot(t,Density,'LineWidth',1.5,'Color',rgb('Yellow'));
    elseif i==4
        plot(t,Density,'LineWidth',1.5,'Color',rgb('Navy'));
    elseif i==5
        plot(t,Density,'LineWidth',1.5,'Color',rgb('Violet'));
    else
        plot(t,Density,'LineWidth',1.5,'Color',rgb('Red'));
    end
    if i==1
      h2=figure('Name','recon rate');
    else
      figure(h2)
      hold on
    end
    %semilogy(t,recon_rate,'LineWidth',1.5)
    %plot(t,recon_rate)
    if i==1
        semilogy(t,recon_rate,'LineWidth',1.5,'Color',rgb('LimeGreen'));
    elseif i==2
        semilogy(t,recon_rate,'LineWidth',1.5,'Color',rgb('Orange'));
    elseif i==3
        semilogy(t,recon_rate,'LineWidth',1.5,'Color',rgb('Yellow'));
    elseif i==4
        semilogy(t,recon_rate,'LineWidth',1.5,'Color',rgb('Navy'));
    elseif i==5
        semilogy(t,recon_rate,'LineWidth',1.5,'Color',rgb('Violet'));
    else
        semilogy(t,recon_rate,'LineWidth',1.5,'Color',rgb('Red'));
    end
    if i==1
      h2b=figure('Name','scaled recon rate');
    else
      figure(h2b)
      hold on
    end
    %semilogy(t,recon_rate,'LineWidth',1.5)
    %plot(t,recon_rate)
    if i==1
        semilogy(t,recon_rate./L,'LineWidth',1.5,'Color',rgb('LimeGreen'));
    elseif i==2
        semilogy(t,recon_rate./L,'LineWidth',1.5,'Color',rgb('Orange'));
    elseif i==3
        semilogy(t,recon_rate./L,'LineWidth',1.5,'Color',rgb('Yellow'));
    elseif i==4
        semilogy(t,recon_rate./L,'LineWidth',1.5,'Color',rgb('Navy'));
    elseif i==5
        semilogy(t,recon_rate./L,'LineWidth',1.5,'Color',rgb('Violet'));
    else
        semilogy(t,recon_rate./L,'LineWidth',1.5,'Color',rgb('Red'));
    end
    if i==1
      h3=figure('Name','mean curvature');
    else
      figure(h3)
      hold on
    end
    %plot(t,mean_curv,'LineWidth',1.5);
    if i==1
        plot(t,mean_curv,'LineWidth',1.5,'Color',rgb('LimeGreen'));
    elseif i==2
        plot(t,mean_curv,'LineWidth',1.5,'Color',rgb('Orange'));
    elseif i==3
        plot(t,mean_curv,'LineWidth',1.5,'Color',rgb('Yellow'));
    elseif i==4
        plot(t,mean_curv,'LineWidth',1.5,'Color',rgb('Navy'));
    elseif i==5
        plot(t,mean_curv,'LineWidth',1.5,'Color',rgb('Violet'));
    else
        plot(t,mean_curv,'LineWidth',1.5,'Color',rgb('Red'));
    end
    if i==1
      h4=figure('Name','LxLyLz');
    else
      figure(h4)
      hold on
    end
    %plot(t_L,total_Lx./total_Lz,'LineWidth',1.5);
    if i==1
        plot(t_L,total_Lx./total_Lz,'LineWidth',1.5,'Color',rgb('LimeGreen'));
    elseif i==2
        plot(t_L,total_Lx./total_Lz,'LineWidth',1.5,'Color',rgb('Orange'));
    elseif i==3
        plot(t_L,total_Lx./total_Lz,'LineWidth',1.5,'Color',rgb('Yellow'));
    elseif i==4
        plot(t_L,total_Lx./total_Lz,'LineWidth',1.5,'Color',rgb('Navy'));
    elseif i==5
        plot(t_L,total_Lx./total_Lz,'LineWidth',1.5,'Color',rgb('Violet'));
    else
        plot(t_L,total_Lx./total_Lz,'LineWidth',1.5,'Color',rgb('Red'));
    end
    if i==1
      h5=figure('Name','Anisotropy');
    else
      figure(h5)
      hold on
    end
    if i==1
        plot(t_Aniso,lpara./lperp,'LineWidth',1.5,'Color',rgb('LimeGreen'));
    elseif i==2
        plot(t_Aniso,lpara./lperp,'LineWidth',1.5,'Color',rgb('Orange'));
    elseif i==3
        plot(t_Aniso,lpara./lperp,'LineWidth',1.5,'Color',rgb('Yellow'));
    elseif i==4
        plot(t_Aniso,lpara./lperp,'LineWidth',1.5,'Color',rgb('Navy'));
    elseif i==5
        plot(t_Aniso,lpara./lperp,'LineWidth',1.5,'Color',rgb('Violet'));
    else
        plot(t_Aniso,lpara./lperp,'LineWidth',1.5,'Color',rgb('Red'));
    end
    clear b
    rr_scaling(counter)=mean(recon_rate(floor(0.1*length(recon_rate)):end))
    D_scaling(counter)=mean(Density(floor(0.1*length(recon_rate)):end))
    counter=counter+1;
end
figure(h1)
xlabel('t','FontSize',16)
ylabel('line density','FontSize',16)
set(gca,'FontSize',16)
hold off
figure(h2)
xlabel('t','FontSize',16)
ylabel('recon rate','FontSize',16)
set(gca,'FontSize',16)
hold off
figure(h2b)
xlabel('t','FontSize',16)
ylabel('recon rate/\Lambda','FontSize',16)
set(gca,'FontSize',16)
hold off
figure(h3)
xlabel('t','FontSize',16)
ylabel('mean curv','FontSize',16)
set(gca,'FontSize',16)
hold off
figure(h4)
xlabel('t','FontSize',16)
ylabel('Lx/Lz','FontSize',16)
set(gca,'FontSize',16)
hold off
figure(h5)
xlabel('t','FontSize',16)
ylabel('lpara/lperp','FontSize',16)
set(gca,'FontSize',16)
hold off
figure
plot(log(D_scaling),log(rr_scaling),'o')
