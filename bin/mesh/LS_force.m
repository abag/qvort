function LS_force(mesh)
close all
load data/dims.log;
load data/LS_forcing.log;
k=LS_forcing(:,1:3);
A=LS_forcing(:,4:6);
B=LS_forcing(:,7:9);
figure('Name','k,A,B')
subplot(1,3,1) ; plot3(k(:,1),k(:,2),k(:,3),'*')
subplot(1,3,2) ; plot3(A(:,1),A(:,2),A(:,3),'*')
subplot(1,3,3) ; plot3(B(:,1),B(:,2),B(:,3),'*')
%define mesh
bsize=dims(2);
for i=1:mesh
  x(i)=(2*i-1)/(2*mesh)*bsize-bsize/2.;
end
vel(mesh,mesh,3)=0.;
for l=1:length(k)
  normk=sqrt(k(l,1)^2+k(l,2)^2+k(l,3)^2);
  for i=1:mesh
    for j=1:mesh
      vel(j,i,:)=cross(A(l,:),k(l,:))*sin(dot(k(l,:),[x(i) x(j) 0.]))+cross(B(l,:),k(l,:))*sin(dot(k(l,:),[x(i) x(j) 0.]));
    end
  end
  vel=vel/normk;
  vel2=sqrt(vel(:,:,1).^2+vel(:,:,2).^2+vel(:,:,3).^2);
  close all
  figure('visible', 'off')
  pcolor(vel2) ; shading interp
  colorbar
  hold on
  quiver(vel(:,:,1),vel(:,:,2),'w')
  print('-dpng',sprintf('data/LS_force%03d',l)) 
end


