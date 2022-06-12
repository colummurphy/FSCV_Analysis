%Break up the groups based on TRT (Target reation time) quartiles for
%each condition (e.g. after break or after big reward)
%Quickly plot rasters and mean traces for each session according to
%conditions defined in call to tr_raster

sessnum=92; %CHANGE HERE
%analyzeSes(sessnum,'fscvchs',[2 4],'nlxch',{'eyex','eyed','pulse','lickx'})
numtargetbreaktrials=length(find(contains({trlists.trlist.type},'targetbreak')));
numalltrials=length(trlists.trlist);
numfixbreaktrials=length(find(contains({trlists.trlist.type},'fixbreak')));
percenttargetbreak=numtargetbreaktrials/numalltrials*100
percentfixbreak=numfixbreaktrials/numalltrials*100
tridsT=[];
tridsF=[];
tridsR=[];
for chnum=[2 4] %CHANGE HERE
%chnum=2;%Change here to look at different sites for given session
sitename=trlists.fscvsites(find([trlists.fscvsites.ch]==chnum)).site;%Get sitename from list
time=[-12,4];%time to extract relative to alignment event, .e.g display_fix
fs=10;%samp freq
tplot=time(1):1/10:time(2); %timestamps in s to plot
alignSamp=find(tplot==0);   %Alignment time point, in samples
%Raster plot after error, target break trials
[axa,tridsT,dataT]=tr_raster(trlists,'fscv','da',chnum,'ttypes',{{'big','left'},{'post','break'}},'win',time,'event','display_fix','sort','rt_target');
close(gcf)
%Raster plot after error, fix break trials
%[axa,tridsF,dataF]=tr_raster(trlists,'fscv','da',chnum,'ttypes',{{'big','left'},{'post','fixbreak'}},'win',time,'event','display_fix','sort','rt_target');
%close(gcf)

%Raster plot after reward ("big" or "small") trials
[axa,tridsR,dataR]=tr_raster(trlists,'fscv','da',chnum,'ttypes',{{'big','left'},{'post','big'}},'win',time,'event','display_fix','sort','rt_target');
close(gcf)

%%
%Plot mean quartiles
qtr1=round(length(tridsT)*.25);
datameanqtr1=nanmean(dataT(1:qtr1,:),1);
datameanqtr4=nanmean(dataT(end-qtr1:end,:),1);
datase1=nanstd(dataT(1:qtr1,:),0,1)./sqrt(qtr1);
datase4=nanstd(dataT(end-qtr1:end,:),0,1)./sqrt(qtr1);


fmeans=figure('color',[1 1 1]); 
set(0,'CurrentFigure',fmeans);    %set figure handle to current figure
axmeans=gca;
plot(axmeans,tplot,datameanqtr1,'-b')
hold on;
plot(axmeans,tplot,datameanqtr1+datase1,'-k')
plot(axmeans,tplot,datameanqtr1-datase1,'-k')
plot(axmeans,tplot,datameanqtr4,'-g')
plot(axmeans,tplot,datameanqtr4-datase4,'-k')
plot(axmeans,tplot,datameanqtr4+datase4,'-k')
ylims=get(axmeans,'ylim');
hold on; plot(axmeans,[0 0],[ylims],'--k')
text1=text(axmeans,0 ,0,'C','color',[0 0 0]);
line2=plot(axmeans,[1.8 1.8],[ylims],'--k'); %target display
text2=text(axmeans,1.8,0,'PT','color',[0 0 0]);
line2.Color(4)=0.5;
xlabel('Time (s)');
ylabel('Delta DA (nM)');
title(['Post Break, Session ' num2str(sessnum) ' | Site ' sitename])
legend(axmeans,'1st qrt TRT','','','4th qrt TRT','','')


qtr=round(length(tridsR)*.25);
datameanqtr1=nanmean(dataR(1:qtr,:),1);
datameanqtr4=nanmean(dataR(end-qtr:end,:),1);
datase1=nanstd(dataR(1:qtr,:),0,1)./sqrt(qtr);
datase4=nanstd(dataR(end-qtr:end,:),0,1)./sqrt(qtr);

fmeans=figure('color',[1 1 1]); 
set(0,'CurrentFigure',fmeans);    %set figure handle to current figure
axmeans=gca;
plot(axmeans,tplot,datameanqtr1,'-b')
hold on;
plot(axmeans,tplot,datameanqtr1+datase1,'-k')
plot(axmeans,tplot,datameanqtr1-datase1,'-k')
plot(axmeans,tplot,datameanqtr4,'-g')
plot(axmeans,tplot,datameanqtr4-datase4,'-k')
plot(axmeans,tplot,datameanqtr4+datase4,'-k')
ylims=get(axmeans,'ylim');
hold on; plot(axmeans,[0 0],[ylims],'--k')
text1=text(axmeans,0 ,0,'C','color',[0 0 0]);
line2=plot(axmeans,[1.8 1.8],[ylims],'--k'); %target display
text2=text(axmeans,1.8,0,'PT','color',[0 0 0]);
line2.Color(4)=0.5;
xlabel('Time (s)');
ylabel('Delta DA (nM)');
title(['Post Big R, Session ' num2str(sessnum) ' | Site ' sitename])
legend(axmeans,'1st qrt TRT','','','4th qrt TRT','','')


end

figure;
trtsT=[trlists.trlist(tridsT).trt];
histogram(trtsT)
trtsR=[trlists.trlist(tridsR).trt];
hold on
histogram(trtsR)
legend('break','rewarded')
