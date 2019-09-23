function vortex_multiplot(filenumbers,varargin)
figure('visible','off');
for i=filenumbers
    subplot(2,2,1);title('tangle 1')
    vortex_plot(i,varargin{:})
    subplot(2,2,2);title('xz-plane')
    vortex_plot(i,varargin{:});view(90,0)
    subplot(2,2,3);title('xy-plane')
    vortex_plot(i,varargin{:});view(0,90)
    subplot(2,2,4);title('tangle 2')
    vortex_plot(i,varargin{:});view(90,90)
    fOUT=sprintf('data/var%04d.jpeg',i)
    print('-djpeg',fOUT)
end
figure('visible','on');
close all
