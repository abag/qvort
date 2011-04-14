%read in the ts file and plot the line length in both normal and logscale
%takes an input that is the number of sections to divide the plot into and give
%a growth/decay rate
function line_length(fit)
if nargin<1
  fit=-5/3
end
A=load('data/ts.log');
t=A(:,2) ; l=A(:,6) ; 
fitted_line=t.^(fit) ;
factor=sum(l(200:length(l)))/sum(fitted_line(200:length(l)))
fitted_line=fitted_line*factor
loglog(t,l,t,fitted_line)
