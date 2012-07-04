%calculates the power density spectra of line length time series
%removes first an final 5% of time series to avoid transient growth
function omega_spectra
global Wrms LL_inside LL_outside LL_total
A=load('./data/ts.log');
n=length(A) ;
t=A(floor(0.08*n):floor(0.9*n),2) ; l=A(floor(0.08*n):floor(0.9*n),6) ;
W=Wrms(floor(0.08*n):floor(0.95*n)) ;
N=length(l) ;
p = abs(fft(l))/(N/2); %% absolute value of the fft
[p_welch f_welch] = pwelch(l,0.2*N,[],[],1/(t(2)-t(1))); %% absolute value of the fft
[p_welch2 f_welch2] = pwelch(LL_inside(floor(0.08*n):floor(0.9*n),9),0.3*N,[],[],1/(t(2)-t(1))); %% absolute value of the fft
dummy_spect=0.2*f_welch.^(-5/3);
loglog(f_welch(4:end),p_welch(4:end),f_welch(6:end),dummy_spect(6:end),f_welch2(3:end),p_welch2(3:end));
return
p = p(1:floor(N/2)).^2;
pW = abs(fft(W))/(N/2); %% absolute value of the fft
pW = pW(1:floor(N/2)).^2;
freq = [0:N/2-1]/t(floor(length(t)));
for i=1:10
  dum_pL=abs(fft(LL_inside(floor(0.08*n):floor(0.9*n),i)))/(N/2);
  pL(i,:)=dum_pL(1:floor(N/2)).^2;
end
%%%%%%%%%%%%%%%%%%%%FITTING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
freq2=freq(floor(0.2*length(freq)):length(freq));
p2=p(floor(0.2*length(p)):length(p));
polyfit(log(freq2),log(p2'),1);
disp(sprintf('spectral slope: %f',ans(1)))
%dummy_spect=freq.^(ans(1));
dummy_spect=freq.^(-5/3);
dummy_spect2=freq.^(1/3);
scaling_factor=p(floor(0.2*length(freq)))/dummy_spect(floor(0.2*length(freq)));
dummy_spect=dummy_spect*scaling_factor;
%%%%%%%%%%%%%%%%%%%%PLOTTING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  figure
  loglog(freq,p,'-m','LineWidth',1.5)
  hold on
  loglog(freq,pW,'-g','LineWidth',1.5)
  set(gca,'FontSize',14)
  xlabel('log freq','FontSize',14)
  ylabel('log power','FontSize',14)
  hold on
  loglog(freq,dummy_spect,'--k','LineWidth',2)
  figure
  cmap=jet(10);
  loglog(freq(2:end-50),smooth(p(2:end-50),2.9),'-','LineWidth',1.5,'Color','b')
  hold on
  loglog(freq(2:end-50),smooth(pL(1,(2:end-50))/100,1),'-','LineWidth',1.5,'Color','r')
  loglog(freq,dummy_spect,'--k','LineWidth',2)
  loglog(freq,dummy_spect2/1E7,'--k','LineWidth',2)

  
