function random_loops_hist
figure('Name','Initial loop distribution')
A=load('data/random_loop_sizes.log');
subplot(2,1,1)
hist(A(:,1))
xlabel('pcount','FontSize',14)
ylabel('N','FontSize',14)
set(gca,'FontSize',14)
subplot(2,1,2)
[bandwidth, n, xout]=kde(A(:,2),2^3) ;
plot(xout,n,'k','LineWidth',2)
axis tight
xlabel('raidus','FontSize',14)
ylabel('PDF','FontSize',14)
set(gca,'FontSize',14)