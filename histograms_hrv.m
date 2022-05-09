event='targ';
trtypes={'big','small'};
targsig='pulse';
sessid='114';
sigtrace={};
mididx=0;
fs=10;
win=[-1 11];  % in seconds starting from aln point
win=[0 4];  % in seconds starting from aln point

hr={};
rmssd={};
rr={};
binwidth=.5;              %ms
binmax=15;
binmin=0;
binedges=binmin:binwidth:binmax;
binwidthrrdev=1;              %ms
binmaxrrdev=27;
binminrrdev=2;
binedgesrrdev=binminrrdev:binwidthrrdev:binmaxrrdev;
binwidthrr=10;              %ms
binmaxrr=600;
binminrr=350;
binedgesrr=binminrr:binwidthrr:binmaxrr;

binwidthhr=2.5;              %ms
binmaxhr=160;
binminhr=100;

binedgeshr=binminhr:binwidthhr:binmaxhr;

binhr={};
binrrstd={};
hrmean={};
binrr={};
for ix=1:length(trtypes)
  targrowt=find(contains({xbinfos.sitelfp},targsig) & ...
      strcmp({xbinfos.event},event) & ...
      contains({xbinfos.sessiontype},trtypes{ix}) & ...
  contains({xbinfos.sessionid},sessid));
    sigtrace{ix}=getfield(xbinfos(targrowt(1)).daall,'lfptracesaln');
    mididx=getfield(xbinfos(targrowt(1)).daall,'mididx');
    winids=[win(1):1/fs:win(2)].*fs+mididx;
    hrtemp=sigtrace{ix};
    rrtemp=1./hrtemp.*60;
    %rmssdtemp=rms(rrtemp);
    hr{ix}=hrtemp(:,winids);
    rr{ix}=rrtemp(:,winids).*1e3;
    rmssd{ix}=rms(diff(rr{ix},1,2),2);
    [binhrv{ix},bb]=histc(rmssd{ix},binedges);
    rrmean{ix}=nanmean(rr{ix},2);
    rrstd{ix}=nanstd(rr{ix},[],2);
    hrmean{ix}=nanmean(hr{ix},2);
    [binrrstd{ix},bb]=histc(rrstd{ix},binedgesrrdev);
    [binhr{ix},bb]=histc(hrmean{ix},binedgeshr);
    [binrr{ix},bb]=histc(rrmean{ix},binedgesrr);
end
  
    fighist=figure('visible','on');
    set(0,'CurrentFigure',fighist); 
    xlims=[binmin binmax];
    subhist(1)=subplot(4,2,1);
    set(fighist, 'Color', [1 1 1]);
    set(fighist,'Position',[50,50,1200,700]);
    countsub=1;
    subhist(1)=subplot(2,2,1);
    subhist(2)=subplot(2,2,2);
        subhist(3)=subplot(2,2,3);
        subhist(4)=subplot(2,2,4);

   bbc= bar(subhist(countsub),binedges,binhrv{1},'FaceColor',[0 0 0],'BarWidth',1,'FaceAlpha',.4,'linestyle','none');
   hold(subhist(countsub),'on');
      bbs= bar(subhist(countsub),binedges,binhrv{2},'FaceColor',[1 0 0],'BarWidth',1,'FaceAlpha',.4,'linestyle','none');
legend(subhist(countsub),'big', 'small','box','off');
    xlabel(subhist(countsub),'RMSSD (ms)');
    set(subhist(countsub),'xtick',binedges(1:10:end));
    set(subhist(countsub),'xlim',xlims,'box','off');
    ylabel(subhist(countsub), '# trials counted')
    set(findall(subhist(countsub),'-property','FontSize'),'FontSize',14)

    countsub=countsub+1;
        
   bbc= bar(subhist(countsub),binedgesrrdev,binrrstd{1},'FaceColor',[0 0 0],'BarWidth',1,'FaceAlpha',.4,'linestyle','none');
   hold(subhist(countsub),'on');
   bbs= bar(subhist(countsub),binedgesrrdev,binrrstd{2},'FaceColor',[1 0 0],'BarWidth',1,'FaceAlpha',.4,'linestyle','none');
legend(subhist(countsub),'big', 'small','box','off');
    set(subhist(countsub),'xtick',binedgesrrdev(1:10:end));
    set(subhist(countsub),'xlim',[binedgesrrdev(1) binedgesrrdev(end)],'box','off');
    xlabel(subhist(countsub),'rr std (ms)');
    ylabel(subhist(countsub), '# trials counted')
        set(findall(subhist(countsub),'-property','FontSize'),'FontSize',14)        
    countsub=countsub+1;
    
    bbc= bar(subhist(countsub),binedgeshr,binhr{1},'FaceColor',[0 0 0],'BarWidth',1,'FaceAlpha',.4,'linestyle','none');
   hold(subhist(countsub),'on');
      bbs= bar(subhist(countsub),binedgeshr,binhr{2},'FaceColor',[1 0 0],'BarWidth',1,'FaceAlpha',.4,'linestyle','none');
legend(subhist(countsub),'big', 'small','box','off');
    set(subhist(countsub),'xtick',binedgeshr(1:10:end));
  set(subhist(countsub),'xlim',[binedgeshr(1) binedgeshr(end)],'box','off');
    xlabel(subhist(countsub),'hr (bpm)');
    ylabel(subhist(countsub), '# trials counted')
        set(findall(subhist(countsub),'-property','FontSize'),'FontSize',14) 
    countsub=countsub+1;
    
    bbc= bar(subhist(countsub),binedgesrr,binrr{1},'FaceColor',[0 0 0],'BarWidth',1,'FaceAlpha',.4,'linestyle','none');
   hold(subhist(countsub),'on');
      bbs= bar(subhist(countsub),binedgesrr,binrr{2},'FaceColor',[1 0 0],'BarWidth',1,'FaceAlpha',.4,'linestyle','none');
legend(subhist(countsub),'big', 'small','box','off');
    set(subhist(countsub),'xtick',binedgesrr(1:10:end));
  set(subhist(countsub),'xlim',[binedgesrr(1) binedgesrr(end)],'box','off');
    xlabel(subhist(countsub),'rr mean (ms)');
    ylabel(subhist(countsub), '# trials counted')
        set(findall(subhist(countsub),'-property','FontSize'),'FontSize',14)    
                    title(subhist(1),[sessid ' | ' event ' | win ' num2str(win(1)) ' - ' num2str(win(2)) ' s'])

        savename=['hist_pulse_' sessid '_' event 'win_' num2str(win(1)) '_' num2str(win(2)) 's'];
    savefig(fighist,[plotparam.savepath savename]);
    %saveas(figsigrep3,[PathName 'tracesnlxbehav.eps'],'epsc')
    saveas(fighist,[plotparam.savepath savename],'tif')