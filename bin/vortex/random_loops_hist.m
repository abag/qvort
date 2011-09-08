%if given the option print will print to .eps file rather than screen
function random_loops_hist(option)
if nargin==0     
  option='empty';
end
switch option
case 'print'
    disp('will not print to screen but instead to .eps files')
case 'empty'
    otherwise
    disp('incorrect option, aborting script and printing help:')
    help curvature
    return
end
switch option
  case 'print'
    figure('visible','off');
  otherwise
    figure('Name','Initial loop distribution')
end
A=load('./data/random_loop_sizes.log');
subplot(2,1,1)
hist(A(:,1))
xlabel('pcount','FontSize',14)
ylabel('N','FontSize',14)
set(gca,'FontSize',14)
subplot(2,1,2)
%[bandwidth, n, xout]=kde(A(:,2),2^3) ;
[n xout]=ksdensity(A(:,2));
plot(xout,n,'k','LineWidth',2)
axis tight
xlabel('radius','FontSize',14)
ylabel('PDF','FontSize',14)
set(gca,'FontSize',14)
switch option
  case 'print'
    disp('printing to random_loop_hist.eps')
    print('-depsc','./random_loop_hist.eps')
end
