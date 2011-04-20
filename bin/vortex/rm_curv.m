%plot a PDF of the curvature of points that have been removed from
%the filament, useful if you have set phone_emission T in run.in
function rm_curv
A=load('./data/removed_curv.log');
B=load('data/curvature.log');
cmax=B(:,3) ; cmax_bar=mean(cmax) ; cmax_latest=cmax(length(cmax)) ;
dims=load('./data/dims.log');
delta=dims(1);
[n xout]=histnorm(A);
dummy_n=min(n):0.000001:max(n);
dummy_xout1(1:length(dummy_n))=sqrt(3)/delta;
dummy_xout2(1:length(dummy_n))=cmax_bar;
dummy_xout3(1:length(dummy_n))=cmax_latest;
plot(xout,n,'LineWidth',2,'Color',rgb('DarkRed'));
hold on
plot(dummy_xout1,dummy_n,'--k',dummy_xout2,dummy_n,'--b',dummy_xout3,dummy_n,'--g');
set(gca,'FontSize',14)
xlabel('\kappa','FontSize',14)
ylabel('PDF(\kappa)','FontSize',14)
