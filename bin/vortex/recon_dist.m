function recon_dist(filenumber)
filename=sprintf('./data/recon_dist%04d.log',filenumber);
fid=fopen(filename);
if fid<0
  disp('file does not exist, exiting script')
  return
end
%read the time
tline=fgetl(fid);
dummy=textscan(tline, '%f');
time=dummy{:}
tline=fgetl(fid);
dummy=textscan(tline, '%f');
angle=dummy{:}
%set ndist
ndist=200;
for j=1:ndist
  tline=fgetl(fid);
  dummy=textscan(tline, '%f');
  dummy_vect=dummy{:};
  dist(j)=dummy_vect(1,:);
  curvi(j)=dummy_vect(2,:);
  curvj(j)=dummy_vect(3,:);
end
plot(dist.^2)