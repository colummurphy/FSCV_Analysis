  trialname=trialtypes.names{it};
binwidth=25;              %nm
binmax=600;
binedges=0:binwidth:binmax;
binsrts=zeros(1,length(binedges));
binsrtsright=zeros(1,length(binedges));
lrts=binfos(81).target_lrt;
rrts=binfos(81).target_rrt;
frts=binfos(81).fix_rt;
[binsrts,bb]=histc(lrts.*1000,binedges);
[binsrtsr,bb]=histc(rrts.*1000,binedges);
slrts=binfos(84).target_lrt;
srrts=binfos(84).target_rrt;
sfrts=binfos(84).fix_rt;
[sbinsrts,bb]=histc(slrts.*1000,binedges);
[sbinsrtsr,bb]=histc(srrts.*1000,binedges);
    fighist=figure('visible','on');
    set(0,'CurrentFigure',fighist); 
    xlims=[0 binmax];
    subhist(1)=subplot(4,2,1);
    set(fighist, 'Color', [1 1 1]);
    set(fighist,'Position',[50,50,1200,500]);
    countsub=1;
    subhist(1)=subplot(1,2,1);
    subhist(2)=subplot(1,2,2);
   bbc= bar(subhist(countsub),binedges,binsrts,'FaceColor',[0 0 0],'BarWidth',1,'FaceAlpha',.4,'linestyle','none');
   hold(subhist(countsub),'on');
      bbs= bar(subhist(countsub),binedges,sbinsrts,'FaceColor',[1 0 0],'BarWidth',1,'FaceAlpha',.4,'linestyle','none');
legend(subhist(countsub),'big', 'small','box','off');
    xlabel(subhist(countsub),'reaction time (left) (ms)');
    set(subhist(countsub),'xtick',binedges(1:10:end));
    set(subhist(countsub),'xlim',xlims,'box','off');
    ylabel(subhist(countsub), '# trials counted')
    countsub=countsub+1;
        set(findall(subhist(countsub),'-property','FontSize'),'FontSize',14)

    bbc= bar(subhist(countsub),binedges,binsrtsr,'FaceColor',[0 0 0],'BarWidth',1,'FaceAlpha',.4,'linestyle','none');
   hold(subhist(countsub),'on');
      bbs= bar(subhist(countsub),binedges,sbinsrtsr,'FaceColor',[1 0 0],'BarWidth',1,'FaceAlpha',.4,'linestyle','none');
legend(subhist(countsub),'big', 'small','box','off');
    xlabel(subhist(countsub),'reaction time (right) (ms)');
    set(subhist(countsub),'xtick',binedges(1:10:end));
    set(subhist(countsub),'xlim',xlims,'box','off');
    ylabel(subhist(countsub), '# trials counted')
        set(findall(subhist(countsub),'-property','FontSize'),'FontSize',14)

     [fbinsrts,bb]=histc(frts.*1000,binedges);
[fbinsrtsr,bb]=histc(sfrts.*1000,binedges);
    bbc= bar(subhist(countsub),binedges,fbinsrts,'FaceColor',[0 0 0],'BarWidth',1,'FaceAlpha',.4,'linestyle','none');
   hold(subhist(countsub),'on');
      bbs= bar(subhist(countsub),binedges,fbinsrtsr,'FaceColor',[1 0 0],'BarWidth',1,'FaceAlpha',.4,'linestyle','none');
legend(subhist(countsub),'big', 'small','box','off');
    xlabel(subhist(countsub),'reaction time (ms)');
    set(subhist(countsub),'xtick',binedges(1:10:end));
    set(subhist(countsub),'xlim',xlims,'box','off');
    ylabel(subhist(countsub), '# trials counted')
        set(findall(subhist(countsub),'-property','FontSize'),'FontSize',14)
        
        
    
    savefig(fighist,[pathrts trialname]);
    %saveas(figsigrep3,[PathName 'tracesnlxbehav.eps'],'epsc')
    saveas(fighist,[pathrts trialname],'tif')