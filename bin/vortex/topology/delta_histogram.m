function delta_histogram(filenumber)
filename=sprintf('./data/delta_adapt%04d.dat',filenumber);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid=fopen(filename);
if fid<0
  disp('var file does not exist, exiting script')
  return
end
time=fread(fid,1,'float64');
number_of_particles=fread(fid,1,'int');
delta=fread(fid,number_of_particles,'float64');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[f,xi] = histnorm(delta);
xlabel('\delta','FontSize',14);
ylabel('PDF(\delta)','FontSize',14);
plot(xi,f,'LineWidth',2)
set(gca,'FontSize',14)