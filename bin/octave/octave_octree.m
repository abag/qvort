function octave_octree(filenumber)
filename=sprintf('data/tree_mesh%03d.log',filenumber);
g=load(filename);
s=size(g);
xmin=min(g(:,1));
xmax=max(g(:,2));
ymin=min(g(:,3));
ymax=max(g(:,4));
zmin=min(g(:,5));
zmax=max(g(:,6));
for i=1:s(1)
  vert1(1)=g(i,1) ; vert1(2)=g(i,3) ; vert1(3)=g(i,5) ;
  vert2(1)=g(i,1) ; vert2(2)=g(i,4) ; vert2(3)=g(i,5) ;
  vert3(1)=g(i,2) ; vert3(2)=g(i,4) ; vert3(3)=g(i,5) ;
  vert4(1)=g(i,2) ; vert4(2)=g(i,3) ; vert4(3)=g(i,5) ;
  vert5(1)=g(i,1) ; vert5(2)=g(i,3) ; vert5(3)=g(i,6) ;
  vert6(1)=g(i,1) ; vert6(2)=g(i,4) ; vert6(3)=g(i,6) ;
  vert7(1)=g(i,2) ; vert7(2)=g(i,4) ; vert7(3)=g(i,6) ;
  vert8(1)=g(i,2) ; vert8(2)=g(i,3) ; vert8(3)=g(i,6) ;
  plot3([vert1(1) vert2(1)],[vert1(2) vert2(2)],[vert1(3) vert2(3)],'k-')
  plot3([vert1(1) vert4(1)],[vert1(2) vert4(2)],[vert1(3) vert4(3)],'k-')
  plot3([vert1(1) vert5(1)],[vert1(2) vert5(2)],[vert1(3) vert5(3)],'k-')
  plot3([vert5(1) vert6(1)],[vert5(2) vert6(2)],[vert5(3) vert6(3)],'k-')
  plot3([vert5(1) vert8(1)],[vert5(2) vert8(2)],[vert5(3) vert8(3)],'k-')
  plot3([vert2(1) vert6(1)],[vert2(2) vert6(2)],[vert2(3) vert6(3)],'k-')
  plot3([vert2(1) vert3(1)],[vert2(2) vert3(2)],[vert2(3) vert3(3)],'k-')
  plot3([vert3(1) vert4(1)],[vert3(2) vert4(2)],[vert3(3) vert4(3)],'k-')
  plot3([vert3(1) vert7(1)],[vert3(2) vert7(2)],[vert3(3) vert7(3)],'k-')
  plot3([vert7(1) vert8(1)],[vert7(2) vert8(2)],[vert7(3) vert8(3)],'k-')
  plot3([vert7(1) vert6(1)],[vert7(2) vert6(2)],[vert7(3) vert6(3)],'k-')  
  plot3([vert8(1) vert4(1)],[vert8(2) vert4(2)],[vert8(3) vert4(3)],'k-')
  hold on
end
hold off
axis([xmin xmax ymin ymax zmin zmax]);
