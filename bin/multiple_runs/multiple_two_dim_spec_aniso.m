function two_dim_spec_aniso(runs,percent,cutoff,fit)
if nargin<3
    do_fit=0;
    cutoff=0;
elseif nargin<4
    do_fit=0;
else
    do_fit=1;
end
cmap=colormap(lines(length(runs)));
color_count=1;
for irun=runs
     if irun==4 || irun==5
        continue
    end
    %generate run prefix
    run_prefix=sprintf('./CF_Tree_spec%01d/',irun);
    %now get filenumbers
    ts=load(strcat(run_prefix,'data/ts.log'));
    %this mst be set to match your run.in file
    mesh_shots_factor=10;
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
        energyr_para=real(fux).^2;
        energyi_para=imag(fux).^2;
        energyr_perp=real(fuy).^2+real(fuz).^2;
        energyi_perp=imag(fuy).^2+imag(fuz).^2;
        midpt=n/2+1;
        spect_para(1:1.5*n)=0.;
        spect_perp(1:1.5*n)=0.;
        if ifile==filenumbers(1)
            avg_spect_para(1:1.5*n)=0.;
            avg_spect_perp(1:1.5*n)=0.;
        end
        for i=1:n
            for j=1:n
                ii=i;
                jj=j;
                if ii>midpt ; ii=n-ii+1; ; end ;
                if jj>midpt ; jj=n-jj+1; ; end ;
                r=int16(sqrt(ii^2+jj^2));
                spect_para(r)=spect_para(r)+3*(energyr_para(i,j)+energyi_para(i,j));
                spect_perp(r)=spect_perp(r)+1.5*(energyr_perp(i,j)+energyi_perp(i,j));
                avg_spect_para(r)=avg_spect_para(r)+3*(energyr_para(i,j)+energyi_para(i,j));
                avg_spect_perp(r)=avg_spect_perp(r)+1.5*(energyr_perp(i,j)+energyi_perp(i,j));
            end
        end
        k=(1:midpt)*(2*pi/dims(2));
    end
    %%%%%%%%%%%%AVERAGE SPECTRA%%%%%%%%%%%%%%%
    avg_spect_perp=avg_spect_perp./length(filenumbers);
    avg_spect_para=avg_spect_para./length(filenumbers);
    if irun==runs(1)
      h1=figure('Name','E(k) parallel')
    else
       figure(h1)
    end
    loglog(k(1:midpt-cutoff),spect_para(1:midpt-cutoff),'Color',cmap(color_count,:),'LineWidth',1.5)
    intervortex=mean(1./sqrt(ts(floor(percent*length(ts)):length(ts),6)/dims(2)^3));
    dumx(1:2)=2*pi/intervortex;
    range=ylim;
    dumy(1)=range(1);
    dumy(2)=range(2);
    hold on
    plot(dumx,dumy,'Color',cmap(color_count,:),'LineWidth',1.5)
    if do_fit==1
        disp(sprintf('fitting a slope of %f to average spect',fit))
        dummy_spect=k.^(fit);
        scaling_factor=0.8*sum(avg_spect_para(1:10))/sum(dummy_spect(1:10));
        dummy_spect=dummy_spect*scaling_factor;
        hold on
        loglog(k,dummy_spect,'--r','LineWidth',2)
    end
    xlabel('log k','FontSize',16) ; ylabel('log E(k)','FontSize',16)
    set(gca,'FontSize',16)
    if irun==runs(1)
      h2=figure('Name','E(k) perpendicular')
    else
       figure(h2)
    end
    loglog(k(1:midpt-cutoff),spect_perp(1:midpt-cutoff),'Color',cmap(color_count,:),'LineWidth',1.5)
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
        scaling_factor=0.8*sum(avg_spect_perp(1:10))/sum(dummy_spect(1:10));
        dummy_spect=dummy_spect*scaling_factor;
        hold on
        loglog(k,dummy_spect,'--r','LineWidth',2)
    end
    xlabel('log k','FontSize',16) ; ylabel('log E(k)','FontSize',16)
    set(gca,'FontSize',16)
end
