function octree(filenumber)
filename=sprintf('data/tree_mesh%03d.log',filenumber);
g=load(filename);
s=size(g);
xmin=min(g(:,1));
xmax=max(g(:,2));
ymin=min(g(:,3));
ymax=max(g(:,4));
zmin=min(g(:,5));
zmax=max(g(:,6));
 fac = [1 2 3 4; ... 
    2 6 7 3; ... 
    4 3 7 8; ... 
    1 5 8 4; ... 
    1 2 6 5; ... 
    5 6 7 8];
for i=1:s(1)
  vert(1,1)=g(i,1) ; vert(1,2)=g(i,3) ; vert(1,3)=g(i,5) ;
  vert(2,1)=g(i,1) ; vert(2,2)=g(i,4) ; vert(2,3)=g(i,5) ;
  vert(3,1)=g(i,2) ; vert(3,2)=g(i,4) ; vert(3,3)=g(i,5) ;
  vert(4,1)=g(i,2) ; vert(4,2)=g(i,3) ; vert(4,3)=g(i,5) ;
  vert(5,1)=g(i,1) ; vert(5,2)=g(i,3) ; vert(5,3)=g(i,6) ;
  vert(6,1)=g(i,1) ; vert(6,2)=g(i,4) ; vert(6,3)=g(i,6) ;
  vert(7,1)=g(i,2) ; vert(7,2)=g(i,4) ; vert(7,3)=g(i,6) ;
  vert(8,1)=g(i,2) ; vert(8,2)=g(i,3) ; vert(8,3)=g(i,6) ;
  patch('Faces',fac,'Vertices',vert,'FaceColor','w');  % patch function
  alpha('clear')
  hold on
end
axis([xmin xmax ymin ymax zmin zmax]);
lighting phong 
view(30,30);