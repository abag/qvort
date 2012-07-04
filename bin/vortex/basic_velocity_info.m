A=load('./data/basic_velocity_info.log');
t=A(:,1);
maxu=A(:,2);
meanu=A(:,3);
minu=A(:,4);
meanux=A(:,5);
meanuy=A(:,6);
meanuz=A(:,7);
figure
subplot(3,1,1)
plot(t,minu,'LineWidth',2)
xlabel('t','FontSize',16)
ylabel('min(u)','FontSize',16)
set(gca,'FontSize',16)
subplot(3,1,2)
plot(t,meanu,'LineWidth',2)
xlabel('t','FontSize',16)
ylabel('mean(u)','FontSize',16)
set(gca,'FontSize',16)
subplot(3,1,3)
plot(t,maxu,'LineWidth',2)
xlabel('t','FontSize',16)
ylabel('max(u)','FontSize',16)
set(gca,'FontSize',16)

figure
subplot(3,1,1)
plot(t,meanux,'LineWidth',2)
xlabel('t','FontSize',16)
ylabel('mean(u_x)','FontSize',16)
set(gca,'FontSize',16)
subplot(3,1,2)
plot(t,meanuy,'LineWidth',2)
xlabel('t','FontSize',16)
ylabel('mean(u_y)','FontSize',16)
set(gca,'FontSize',16)
subplot(3,1,3)
plot(t,meanuz,'LineWidth',2)
xlabel('t','FontSize',16)
ylabel('mean(u_z)','FontSize',16)
set(gca,'FontSize',16)
