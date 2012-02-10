function two_dim_spec_aniso(runs,percent,cutoff,fit)
if nargin<3
    do_fit=0;
    cutoff=0;
elseif nargin<4
    do_fit=0;
else
    do_fit=1;
end
cmap=colormap(jet(length(runs)));
color_count=1;
for irun=runs
    %generate run prefix
    run_prefix=sprintf('./run%01d/',irun);
    %now get filenumbers
    ts=load(strcat(run_prefix,'data/ts.log'));
    %this mst be set to match your run.in file
    mesh_shots_factor=20;
    filenumbers=floor(percent*(length(ts)/mesh_shots_factor)):((length(ts)-1)/mesh_shots_factor);
    for ifile=filenumbers
        filename=strcat(run_prefix,sprintf('data/vel_slice_2D%04d.dat',ifile))
        fid=fopen(filename);
        if fid<0
            disp('2D slice file does not exist, exiting script')
            return
        end
        dims=load(strcat(run_prefix,'data/dims.log'));
        A=fread(fid,'float64');
        s=length(A); s=s/8; s=sqrt(s);
        B=reshape(A,8,s,s);
        x=squeeze(B(1,:,:));
        y=squeeze(B(2,:,:));
        xx=squeeze(x(1,1,:));
        yy=squeeze(y(1,:,1));
        usupx=squeeze(B(3,:,:));
        usupy=squeeze(B(4,:,:));
        usupz=squeeze(B(5,:,:));
        unormx=squeeze(B(6,:,:));
        unormy=squeeze(B(7,:,:));
        unormz=squeeze(B(8,:,:));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SPECTRA%%%%%%%%%%%%%%%%%%%
        n=s;
        fux=fftn(usupx)/(n^2);
        fuy=fftn(usupy)/(n^2);
        fuz=fftn(usupz)/(n^2);
        energyr=real(fuy).^2+real(fuz).^2;
        energyi=imag(fuy).^2+imag(fuz).^2;
        midpt=n/2+1;
        spect(1:1.5*n)=0.;
        if ifile==filenumbers(1)
            avg_spect(1:1.5*n)=0.;
        end
        for i=1:n
            for j=1:n
                ii=i;
                jj=j;
                if ii>midpt ; ii=n-ii+1; ; end ;
                if jj>midpt ; jj=n-jj+1; ; end ;
                r=int16(sqrt(ii^2+jj^2));
                spect(r)=spect(r)+energyr(i,j)+energyi(i,j);
                avg_spect(r)=avg_spect(r)+energyr(i,j)+energyi(i,j);
            end
        end
        k=(1:midpt)*(2*pi/dims(2));
    end
    %%%%%%%%%%%%AVERAGE SPECTRA%%%%%%%%%%%%%%%
    avg_spect=avg_spect./length(filenumbers);
    figure('Name','average E(k)')
    loglog(k(1:midpt-cutoff),spect(1:midpt-cutoff),'Color',cmap(color_count,:),'LineWidth',1.5)
    intervortex=mean(1./sqrt(ts(floor(percent*length(ts)):length(ts),6)/dims(2)^3));
    dumx(1:2)=2*pi/intervortex;
    range=ylim;
    dumy(1)=range(1);
    dumy(2)=range(2);
    hold on
    plot(dumx,dumy,'Color',cmap(color_count,:),'LineWidth',1.5)
    color_count=color_count+1;
    if do_fit==1
        disp(sprintf('fitting a slope of %f to average spect',fit))
        dummy_spect=k.^(fit);
        scaling_factor=0.8*sum(avg_spect(1:10))/sum(dummy_spect(1:10));
        dummy_spect=dummy_spect*scaling_factor;
        hold on
        loglog(k,dummy_spect,'--r','LineWidth',2)
    end
    xlabel('log k','FontSize',16) ; ylabel('log E(k)','FontSize',16)
    set(gca,'FontSize',16)
end
