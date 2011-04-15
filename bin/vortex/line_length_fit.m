%read in the ts file and plot the line length in logscale
%takes an input that is a line to fit to log-log plot
function line_length(fit)
if nargin<1
  fit=-5/3
end
A=load('data/ts.log');
t=A(:,2) ; l=A(:,6) ; 
fitted_line=t.^(fit) ;
%get the fitted line to lie over the data
factor=sum(l(floor(0.75*length(l)):length(l)))/sum(fitted_line(floor(0.75*length(l)):length(l)));
fitted_line=1.5*fitted_line*factor;
loglog(t,l,'k*','MarkerSize',4)
hold on
loglog(t,fitted_line,'r-','LineWidth',2)
hold off
set(gca,'FontSize',14)
xlabel('log t','FontSize',14)
ylabel('log L','FontSize',14)
text(sum((t))./size(t),sum((fitted_line))./size(t),num2str(fit),'FontSize',14)
