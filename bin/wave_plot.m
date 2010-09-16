wave=load('./data/wave_info.log');
subplot(1,2,1)
  plot(wave(:,1),wave(:,2),'.-r','LineWidth',2)
  xlabel('k')
  ylabel('Amp')
subplot(1,2,2)
  loglog(wave(:,1),wave(:,2),'-m','LineWidth',2)
  xlabel('log k')
  ylabel('log Amp')
  