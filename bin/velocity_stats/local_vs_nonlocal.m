function local_vs_nonlocal
A=load('./data/local_v_nonlocal.log');
figure('Name','local (induced) contribution')
plot(A(:,1),A(:,3),'g','LineWidth',2)
hold on
plot(A(:,1),A(:,2),'r','LineWidth',2)
hold on
plot(A(:,1),A(:,4),'k','LineWidth',2)
set(gca,'FontSize',16)
xlabel('t','FontSize',16)
ylabel('u_{loc}','FontSize',16)
figure('Name','non-local (BS) contribution')
plot(A(:,1),A(:,6),'g','LineWidth',2)
hold on
plot(A(:,1),A(:,5),'r','LineWidth',2)
hold on
plot(A(:,1),A(:,7),'k','LineWidth',2)
set(gca,'FontSize',16)
xlabel('t','FontSize',16)
ylabel('u_{BS}','FontSize',16)
figure('Name','u_{BS}/u_{LIA}')
plot(A(:,1),A(:,6)./A(:,3),'k','LineWidth',2)
set(gca,'FontSize',16)
xlabel('t','FontSize',16)
ylabel('u_{BS}/u_{LIA}','FontSize',16)
figure('Name','logscale u_{BS}/u_{LIA}')
semilogy(A(:,1),A(:,6)./A(:,3),'k','LineWidth',2)
set(gca,'FontSize',16)
xlabel('t','FontSize',16)
ylabel('u_{BS}/u_{LIA}','FontSize',16)