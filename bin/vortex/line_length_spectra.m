%calculates the power density spectra of line length time series
%removes first an final 5% of time series to avoid transient growth
function line_length_spectrum
A=load('./data/ts.log');
n=length(A) ;
t=A(floor(0.05*n):floor(0.95*n),2) ; l=A(floor(0.05*n):floor(0.95*n),6) ;
N=length(l) ;
p = abs(fft(l))/(N/2); %% absolute value of the fft
p = p(1:floor(N/2)).^2;
freq = [0:N/2-1]/t(floor(length(t)));
%%%%%%%%%%%%%%%%%%%%FITTING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
freq2=freq(floor(0.2*length(freq)):length(freq));
p2=p(floor(0.2*length(p)):length(p));
polyfit(log(freq2),log(p2'),1);
disp(sprintf('spectral slope: %f',ans(1)))
%dummy_spect=freq.^(ans(1));
dummy_spect=freq.^(-5/3);
scaling_factor=p(floor(0.2*length(freq)))/dummy_spect(floor(0.2*length(freq)));
dummy_spect=dummy_spect*scaling_factor;
%%%%%%%%%%%%%%%%%%%%PLOTTING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(3)/1.2 scrsz(4)/2],'PaperPosition',[0.25 2.5 28.0 12.0],'color','w','visible','on','Name', 'line length spectrum')
subplot(1,2,1)
  plot(t,l,'-r','LineWidth',2);
  set(gca,'FontSize',14)
  xlabel('t','FontSize',14)
  ylabel('l','FontSize',14)
subplot(1,2,2)
  loglog(freq,p,'-m','LineWidth',1.5)
  set(gca,'FontSize',14)
  xlabel('log freq','FontSize',14)
  ylabel('log power','FontSize',14)
  hold on
  loglog(freq,dummy_spect,'--k','LineWidth',2)
