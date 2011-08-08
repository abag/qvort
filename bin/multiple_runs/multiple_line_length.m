function multiple_line_length(nruns)
if nargin<1
 nruns=1
end
figure
for i=1:nruns
  filename=sprintf('./run%d/data/ts.log',i);
  filename2=sprintf('./run%d/data/dims.log',i);
  A=load(filename);
  B=load(filename2);
  plot(A(:,2),A(:,6)/B(2)^3)
  hold on
end 
