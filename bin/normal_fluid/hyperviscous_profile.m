A=load('./data/hyperviscous.log');
plot(A(:,1),A(:,2),'-k','LineWidth',2)
xlabel('\kappa','FontSize',16)
ylabel('\alpha','FontSize',16)
set(gca,'FontSize',16)
