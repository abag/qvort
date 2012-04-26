function multiple_recon_fit(start,finish,skip)
global a b angle R ind cmax
for i=start:skip:finish
    [a(i) b(i) c(i) R(i) angle(i) cmax(i)] = recon_dist(i,0);
end
k=sqrt(0.997E-03);
a=a/k;
b=b./(a*k);
figure
subplot(1,2,1)
ind = find(R>0.98 & a>0. & b<10);
ksdensity(a(ind))
subplot(1,2,2)
ksdensity(b(ind))