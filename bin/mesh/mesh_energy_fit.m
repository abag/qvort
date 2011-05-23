function mesh_energy_fit(fit,smoothed)
if nargin<1
  fit=-5/3
  smoothed=1
elseif nargin<2
  smoothed=1
end
load mesh_energy.mat
fitted_line=t.^(fit) ;
%get the fitted line to lie over the data
factor=sum(u2(floor(0.75*length(u2)):length(u2)))/sum(fitted_line(floor(0.75*length(u2)):length(u2)));
fitted_line=1.5*fitted_line*factor;
loglog(t,smooth(u2,smoothed),'k*','MarkerSize',4)
hold on
loglog(t,fitted_line,'r-','LineWidth',2)
hold off
set(gca,'FontSize',14)
xlabel('log t','FontSize',14)
ylabel('log E','FontSize',14)
text(sum((t))./size(t),sum((fitted_line))./size(t),num2str(fit),'FontSize',14)
