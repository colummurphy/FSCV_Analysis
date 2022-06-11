%Quickly plot rasters and mean traces for each session according to
%conditions defined in call to tr_raster
sessnum=92;
%analyzeSes(sessnum,'fscvchs',[2 4],'nlxch',{'eyex','eyed','pulse','lickx'})
chnum=4;%Change here to look at different sites for given session
sitename=trlists.fscvsites(find([trlists.fscvsites.ch]==chnum)).site;%Get sitename from list
time=[-12,4];%time to extract relative to alignment event, .e.g display_fix
fs=10;%samp freq
tplot=time(1):1/10:time(2); %timestamps in s to plot
alignSamp=find(tplot==0);   %Alignment time point, in samples
%Raster plot after error, target break trials
[axa,trids,dataT]=tr_raster(trlists,'fscv','da',chnum,'ttypes',{{'big','left'},{'post','targetbreak'}},'win',time,'event','display_fix');
%Raster plot after error, fix break trials
[axa,trids,dataF]=tr_raster(trlists,'fscv','da',chnum,'ttypes',{{'big','left'},{'post','fixbreak'}},'win',time,'event','display_fix');
%Raster plot after reward ("big" or "small") trials
[axa,trids,dataR]=tr_raster(trlists,'fscv','da',chnum,'ttypes',{{'big','left'},{'post','big'}},'win',time,'event','display_fix');
%%
datamean=nanmean(dataT,1);
datase=nanstd(dataT,0,1)./sqrt(size(dataT,1));
fmeans=figure('color',[1 1 1]); 
set(0,'CurrentFigure',fmeans);    %set figure handle to current figure
axmeans=gca;
plot(axmeans,tplot,datamean,'-b')
hold on;
plot(axmeans,tplot,datamean+datase,'-k')
plot(axmeans,tplot,datamean-datase,'-k')
datamean=nanmean(dataF,1);
datase=nanstd(dataF,0,1)./sqrt(size(dataF,1));
plot(axmeans,tplot,datamean,'-g')
plot(axmeans,tplot,datamean-datase,'-k')
plot(axmeans,tplot,datamean+datase,'-k')

datamean=nanmean(dataR,1);
datase=nanstd(dataR,0,1)./sqrt(size(dataR,1));
plot(axmeans,tplot,datamean,'-r')
plot(axmeans,tplot,datamean-datase,'-k')
plot(axmeans,tplot,datamean+datase,'-k')

ylims=get(axmeans,'ylim');
hold on; plot(axmeans,[0 0],[ylims],'--k')
text1=text(axmeans,0 ,0,'C','color',[0 0 0]);
line2=plot(axmeans,[1.8 1.8],[ylims],'--k'); %target display
text2=text(axmeans,1.8,0,'PT','color',[0 0 0]);
line2.Color(4)=0.5;
xlabel('Time (s)');
ylabel('Delta DA (nM)');
title(['Session ' num2str(sessnum) ' | Site ' sitename])
legend(axmeans,'Post PT Break','','','Post C Break','','','Post Big R')