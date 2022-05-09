savepath=fullfile('Z:\inj-monkey2\analysis\cleopatra\mult', filesep, 'multifft' ,filesep);
if ~isdir(savepath)
    mkdir(savepath)
end
savename=[savepath 'beta_spec_data_'];
cnputfft={};
regs=['cp'];
for ireg=1:2
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
targsites=find(contains({fftdata.sitelfp},regs(ireg)));
targdata=fftdata(targsites);
binedges=10:1:40;
[bindata,bb]=histc([targdata.fmin],binedges);

histplot= bar(axa{1},binedges,bindata,'FaceColor',[0 1 0],'BarWidth',1,'FaceAlpha',.4,'linestyle','none');
hold(axa{1},'on');
[bindata2,bb]=histc([targdata.fmax],binedges);
histplot2= bar(axa{1},binedges,bindata2,'FaceColor',[0 0 1],'BarWidth',1,'FaceAlpha',.25,'linestyle','none');
[bindata3,bb]=histc([targdata.fpeak],binedges);
histplot2= bar(axa{1},binedges,bindata3,'FaceColor',[1 0 0],'BarWidth',1,'FaceAlpha',.35,'linestyle','none');
meanfmin=nanmean([targdata.fmin])
meanfmax=nanmean([targdata.fmax])
meanfpeak=nanmean([targdata.fpeak])
cnputfft(ireg).fmin=meanfmin;
cnputfft(ireg).fmax=meanfmax;
cnputfft(ireg).fpeak=meanfpeak;

ylims=get(axa{1},'ylim');
line1=plot(axa{1},[meanfmin meanfmin],ylims,'color',[0 1 0]);
line1.Color(4)=0.4;
line2=plot(axa{1},[meanfmax meanfmax],ylims,'color',[0 0 1]);
line2.Color(4)=0.4;
line3=plot(axa{1},[meanfpeak meanfpeak],ylims,'color',[1 0 0]);
line3.Color(4)=0.4;
legend(axa{1},'f_L','f_H','f_C','mean f_L','mean f_H', 'mean f_C');
xlabel(axa{1},'Freq. (Hz)');
ylabel(axa{1},'# sites');
if strcmp(regs(ireg),'c')
title(axa{1},['cn session-sites (' num2str(length(targdata)) ')']);
savename=[savename 'cn'];
else
    title(axa{1},['put session-sites (' num2str(length(targdata)) ')']);
    savename=[savename 'put'];

end
cla(axa{2})
binedges=0:1:20;
[bindata_widths,bb]=histc([targdata.fmax]-[targdata.fmin],binedges);
histplot= bar(axa{2},binedges,bindata_widths,'FaceColor',[0 0 0],'BarWidth',1,'FaceAlpha',.4,'linestyle','none');
meanfspan=nanmean([targdata.fmax]-[targdata.fmin])
ylims=get(axa{2},'ylim');
line1=plot(axa{2},[meanfspan meanfspan],ylims,'color',[0 0 0]);
line1.Color(4)=0.4;
xlabel(axa{2},'Beta Freq. Span (Hz)');
ylabel(axa{2},'# sites');
cnputfft(ireg).fspan=meanfspan;

cla(axa{3})
binedges=10:1:40;
targsites=find(contains({avgfft.sitelfp},regs(ireg)));
targdata=avgfft(targsites);
[bindata,bb]=histc([targdata.meanfmin],binedges);
histplot= bar(axa{3},binedges,bindata,'FaceColor',[0 1 0],'BarWidth',1,'FaceAlpha',.4,'linestyle','none');
[bindata2,bb]=histc([targdata.meanfmax],binedges);
hold(axa{3},'on'); 
histplot2= bar(axa{3},binedges,bindata2,'FaceColor',[0 0 1],'BarWidth',1,'FaceAlpha',.25,'linestyle','none');
[bindata3,bb]=histc([targdata.meanfpeak],binedges);
histplot2= bar(axa{3},binedges,bindata3,'FaceColor',[1 0 0],'BarWidth',1,'FaceAlpha',.35,'linestyle','none');
meanfmin=nanmean([targdata.meanfmin])
meanfmax=nanmean([targdata.meanfmax])
meanfpeak=nanmean([targdata.meanfpeak])
ylims=get(axa{3},'ylim');
line1=plot(axa{3},[meanfmin meanfmin],ylims,'color',[0 1 0]);
line1.Color(4)=0.4;
line2=plot(axa{3},[meanfmax meanfmax],ylims,'color',[0 0 1]);
line2.Color(4)=0.4;
line3=plot(axa{3},[meanfpeak meanfpeak],ylims,'color',[1 0 0]);
line3.Color(4)=0.4;
legend(axa{3},'f_L','f_H','f_C','mean f_L','mean f_H', 'mean f_C');
cnputfft(ireg).elecfmin=meanfmin;
cnputfft(ireg).elecfmax=meanfmax;
cnputfft(ireg).elecfpeak=meanfpeak;

%legend(axa{3},'f_L','f_H','f_C');
xlabel(axa{3},'Freq. (Hz)');
ylabel(axa{3},'# sites');
title(axa{3},['unique electrode pairs (' num2str(length(targdata)) ')']);
cla(axa{4})
binedges=0:1:20;
[bindata_widths,bb]=histc([targdata.meanfmax]-[targdata.meanfmin],binedges);
histplot= bar(axa{4},binedges,bindata_widths,'FaceColor',[0 0 0],'BarWidth',1,'FaceAlpha',.4,'linestyle','none');
meanfspan=nanmean([targdata.meanfmax]-[targdata.meanfmin])
ylims=get(axa{4},'ylim');
line1=plot(axa{4},[meanfspan meanfspan],ylims,'color',[0 0 0]);
line1.Color(4)=0.4;
xlabel(axa{4},'Beta Freq. Span (Hz)');
ylabel(axa{4},'# sites');
cnputfft(ireg).elecfspan=meanfspan;



savefig(figsess,savename);
saveas(figsess,savename,'jpg')
print(figsess,savename,'-painters','-depsc');
save(savename,'fftdata','avgfft');
end
