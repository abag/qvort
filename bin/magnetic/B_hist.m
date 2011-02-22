function B_hist(filenumber)
filename=sprintf('data/full_B%04d.dat',filenumber);
fid=fopen(filename);
if fid<0
    disp('var file does not exist, exiting script')
    return
end
number_of_particles=fread(fid,1,'int');
time=fread(fid,1,'float64');
B=fread(fid,number_of_particles,'float64');
B2=B(B>0);
%[bandwidth,density,xmesh]=kde(B2,2^4);
[density,xmesh]=ksdensity(B2);
plot(xmesh,density,'LineWidth',2,'Color','b')
axis tight
xlabel('B','FontSize',16)
ylabel('density','FontSize',16)
set(gca,'FontSize',16)
