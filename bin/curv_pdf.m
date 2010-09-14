%a script to plot histograms of curvature
%note I assume that the number of bins is 10
function curv_hist
A=load('./data/curv_pdf.log');
s=size(A) ; snap_number=s(1)/10 ;
B=reshape(A,10,snap_number,2) ;
cmap=colormap(jet(snap_number)) ;
for i=1:snap_number
 plot(B(:,i,1),B(:,i,2),'-','Color',cmap(i,:)) ;
 hold on   
end
hold off
colorbar