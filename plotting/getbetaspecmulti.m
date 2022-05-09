function [fftdata avgfft]=getbetaspecmulti(plotparam,xinfos,binfos,varargin)
fftdata=[];
win=[-3 4];
freq=[10 50];
markers={'fix','fixeye','targeye','outcome'};
conditions={'big'};
plottype='plotboth';
event='targ';
avgchs=1;
savepath=fullfile(plotparam.savepath, 'multifft' ,filesep);
if ~isdir(savepath)
    mkdir(savepath)
end
savename=[savepath 'beta_spec_data_all'];

argnum=1;
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'win'
            argnum=argnum+1;
            win=varargin{argnum};
        case 'event'
            argnum=argnum+1;
            event=varargin{argnum};     %alignment event
        case 'freq'
            argnum=argnum+1;
            freq=varargin{argnum};
        case 'markers'
            argnum=argnum+1;
            markers=varargin{argnum};
        case 'conditions'
            argnum=argnum+1;
            conditions=varargin{argnum};        %all conditions to include
        case 'avgchs'
            %average same electrode pairs, ie. same lfpchs name
            avgchs=1;
    end
    argnum=argnum+1;
end


%get beta spectrum characteristics, peak, width, low freq, upper freq
%limit, from multiple sessions
%called from plotmultiple.m

validlfps=plotparam.lfpchs;

xids=find(strcmp({xinfos.event},event) & contains({xinfos.sessiontype},'big') & ismember({xinfos.sitelfp},validlfps) & strcmp({xinfos.subj},'patra'));      %still need to remove duplicates, ie 2 da paired with same lfp in xinfos
count=1;
for ix=1:length(xids)
    sitelfp=xinfos(xids(ix)).sitelfp;
    sessnum=str2num(xinfos(xids(ix)).sessionid);
    disp([num2str(ix) '/' num2str(length(xids)) ': sess ' num2str(sessnum) ' sitelfp ' sitelfp]);
    xsessids=[];
    if ~isempty(fftdata)
        xsessids=find([fftdata.sessnum]==sessnum);
    end
    if isempty(fftdata)
        maxfft=plotspectrosel(sessnum,sitelfp,[],event,plotparam,'win',win,'freq',freq,'markers',markers,'xinfos',xinfos,'binfos',binfos,'types',conditions,plottype);
        maxfft.sessnum=sessnum;
        maxfft.sitelfp=sitelfp;
        count=count+1;
        fftdata=[fftdata maxfft];
    elseif isempty(xsessids)
        %new session
        maxfft=plotspectrosel(sessnum,sitelfp,[],event,plotparam,'win',win,'freq',freq,'markers',markers,'xinfos',xinfos,'binfos',binfos,'types',conditions,plottype);
        maxfft.sessnum=sessnum;
        maxfft.sitelfp=sitelfp;
        count=count+1;
        fftdata=[fftdata maxfft];
    elseif ~contains(sitelfp,{fftdata(xsessids).sitelfp})
        %check not duplicate, if already doing session
        maxfft=plotspectrosel(sessnum,sitelfp,[],event,plotparam,'win',win,'freq',freq,'markers',markers,'xinfos',xinfos,'binfos',binfos,'types',conditions,plottype);
        maxfft.sessnum=sessnum;
        maxfft.sitelfp=sitelfp;
        count=count+1;
        fftdata=[fftdata maxfft];
    else
        disp([num2str(ix) '/' num2str(length(xids)) ': sess ' num2str(sessnum) ' sitelfp ' sitelfp ' skipping--duplicate ']);
    end
    close all;
end
avgfft=[];
if avgffts
    uniquelfps=unique({fftdata.sitelfp});
    avgfft=[];
    for il=1:length(uniquelfps)
        curlfp=uniquelfps{il};
        ids=find(contains({fftdata.sitelfp},curlfp));
        fmins=[];
        fmaxs=[];
        fpeaks=[];
        fspans=[];
        fmins=[fftdata(ids).fmin];
        fmaxs=[fftdata(ids).fmax];
        fpeaks=[fftdata(ids).fpeak];
        curfft.sitelfp=curlfp;
        curfft.sessnums=[fftdata(ids).sessnum];
        curfft.fmins=fmins;
        curfft.fmaxs=fmaxs;
        curfft.fpeaks=fpeaks;
        curfft.meanfmin=nanmean(fmins);
        curfft.stdfmin=nanstd(fmins);
        curfft.cifmin=nanstd(fmins)./sqrt(length(fmins))*1.96;
        curfft.meanfpeak=nanmean(fpeaks);
        curfft.stdfpeak=nanstd(fpeaks);
        curfft.cifpeak=nanstd(fpeaks)./sqrt(length(fmins))*1.96;        
        curfft.meanfmax=nanmean(fmaxs);
        curfft.stdfmax=nanstd(fmaxs);
        curfft.cifmax=nanstd(fmaxs)./sqrt(length(fmins))*1.96;        
        avgfft=[avgfft curfft];
    end
    
end

figsess=figure('position',[50 50 1000 800],'color',[1 1 1]);
set(0,'CurrentFigure',figsess);    %set figure handle to current figure
axa={};
clf(figsess,'reset');
set(figsess,'color',[1 1 1]);
axpos={};
for ip=1:4
    axa{ip}=subplot(2,2,ip);   hold(axa{ip},'on');
    set(axa{ip},'units','pixels');
    axpos{ip}=get(axa{ip},'position');    
   % set(axa{ip},'position',[(ip-1)*axsiz(1)+ip*(mar/2)+off,125,axsiz(1),axsiz(2)])
       % ax=text(axa{ip},10,axsiz(2),plotvar{ip},'units','pixels','fontweight','bold');
end
cla(axa{1})
binedges=10:1:40;
[bindata,bb]=histc([fftdata.fmin],binedges);
histplot= bar(axa{1},binedges,bindata,'FaceColor',[0 1 0],'BarWidth',1,'FaceAlpha',.4,'linestyle','none');
hold(axa{1},'on');
[bindata2,bb]=histc([fftdata.fmax],binedges);
histplot2= bar(axa{1},binedges,bindata2,'FaceColor',[0 0 1],'BarWidth',1,'FaceAlpha',.25,'linestyle','none');
[bindata3,bb]=histc([fftdata.fpeak],binedges);
histplot2= bar(axa{1},binedges,bindata3,'FaceColor',[1 0 0],'BarWidth',1,'FaceAlpha',.35,'linestyle','none');
legend(axa{1},'f_L','f_H','f_C');
xlabel(axa{1},'Freq. (Hz)');
ylabel(axa{1},'# sites');
title(axa{1},['all session-sites (' num2str(length(fftdata)) ')']);

cla(axa{2})
binedges=0:1:20;
[bindata_widths,bb]=histc([fftdata.fmax]-[fftdata.fmin],binedges);
histplot= bar(axa{2},binedges,bindata_widths,'FaceColor',[0 0 0],'BarWidth',1,'FaceAlpha',.4,'linestyle','none');
xlabel(axa{2},'Beta Freq. Span (Hz)');
ylabel(axa{2},'# sites');

cla(axa{3})
binedges=10:1:40;
[bindata,bb]=histc([avgfft.meanfmin],binedges);
histplot= bar(axa{3},binedges,bindata,'FaceColor',[0 1 0],'BarWidth',1,'FaceAlpha',.4,'linestyle','none');
[bindata2,bb]=histc([avgfft.meanfmax],binedges);
hold(axa{3},'on'); 
histplot2= bar(axa{3},binedges,bindata2,'FaceColor',[0 0 1],'BarWidth',1,'FaceAlpha',.25,'linestyle','none');
[bindata3,bb]=histc([avgfft.meanfpeak],binedges);
histplot2= bar(axa{3},binedges,bindata3,'FaceColor',[1 0 0],'BarWidth',1,'FaceAlpha',.35,'linestyle','none');
legend(axa{3},'f_L','f_H','f_C');
xlabel(axa{3},'Freq. (Hz)');
ylabel(axa{3},'# sites');
title(axa{3},['unique electrode pairs (' num2str(length(avgfft)) ')']);
cla(axa{4})
binedges=0:1:20;
[bindata_widths,bb]=histc([avgfft.meanfmax]-[avgfft.meanfmin],binedges);
histplot= bar(axa{4},binedges,bindata_widths,'FaceColor',[0 0 0],'BarWidth',1,'FaceAlpha',.4,'linestyle','none');
xlabel(axa{4},'Beta Freq. Span (Hz)');
ylabel(axa{4},'# sites');



savefig(figsess,savename);
saveas(figsess,savename,'jpg')
print(figsess,savename,'-painters','-depsc');
save(savename,'fftdata','avgfft');

