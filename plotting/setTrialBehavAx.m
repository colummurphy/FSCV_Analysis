function rts=setTrialBehavAx(data,plotparam,hfig,seltrials,plotnum)
%place on side of trial by trial data, other behavioral metrics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%behavioral parameters flipped vertical

colorsize=plotparam.colorsize;
vertplotwidth=plotparam.vertplotwidth;
vertplotwidth2=plotparam.vertplotwidth2;
colortrial=plotparam.colormap;
colorm=plotparam.markercolors;
numtrials=size(seltrials,2);
win=plotparam.win;      %in fscv samples
origpos=getpixelposition(hfig{plotnum});
vdata=data.target_rt(seltrials);       %NO invert because plotting downwards, use 'ydir' instead
if strcmp(plotparam.trialtype,'fixbreak')
    vdata=data.fix_rt(seltrials);
end
scaleRT=.05;
plot(hfig{plotnum},vdata,1:numtrials,'ksq','markersize',3)
set(hfig{plotnum},'box','off')
set(hfig{plotnum},'xtick',[],'xcolor',[1 1 1]); pause(0.1);
set(hfig{plotnum},'ytick',[],'ycolor',[1 1 1]);  pause(0.1);

xvert=xlabel(hfig{plotnum},'target rt','color',[0 0 0]); 
if strcmp(plotparam.trialtype,'fixbreak')
    xvert=xlabel(hfig{plotnum},'fix rt','color',[0 0 0]); 
end
set(xvert,'units','pixels');
xvertpos=get(xvert,'position');
%xvertpos=getpixelposition(xvert);
set(xvert,'Units','Pixel','Position', [xvertpos(1)+5 xvertpos(2)+colorsize(2)+50],'rotation',90);
hold(hfig{plotnum},'on')
set(hfig{plotnum},'YDir','reverse') 

if numtrials<=2
    %quit if too little trials
    rts=nan;
    return;
end
%plot scale bar
line([origpos(1)+colorsize(1)+10; origpos(1)+colorsize(1)+10],[origpos(2); origpos(2)])
plot(hfig{plotnum},[.1; .1+scaleRT], [1; 1], '-k','linewidth',1)
text(hfig{plotnum},scaleRT+.04,2.5,'.05 s','color',[0 0 0]);
ylim(hfig{plotnum},[1 numtrials+1]);
xlim(hfig{plotnum},[.1 max(vdata)]);

hold(hfig{plotnum},'off')

%Current trial reward side
vdata2=data.rewardside(seltrials);
plot(hfig{plotnum+1},1,1:numtrials,'color',[0 0 0]);
set(hfig{plotnum+1},'box','off')
set(hfig{plotnum+1},'xtick',[],'xcolor',[1 1 1]); pause(0.1);
set(hfig{plotnum+1},'ytick',[],'ycolor',[1 1 1]);  pause(0.1);
xvert=xlabel(hfig{plotnum+1},'reward side','color',[0 0 0]); 
set(xvert,'units','pixels');
xvertpos=get(xvert,'position');

%xvertpos=getpixelposition(xvert);
set(xvert,'Units','Pixel','Position', [xvertpos(1)+5 xvertpos(2)+colorsize(2)+50],'rotation',90);
hold(hfig{plotnum+1},'on')
for ii=1:numtrials
    if vdata2(ii)==0
        text(hfig{plotnum+1},1,ii,'\leftarrow','fontsize',8)
    else
        text(hfig{plotnum+1},1,ii,'\rightarrow','fontsize',8)
    end
end
ylim(hfig{plotnum+1},[1 numtrials+1])
set(hfig{plotnum+1},'YDir','reverse') 
hold(hfig{plotnum+1},'off')

%Current big reward side
vdata2b=data.bigside(seltrials);
plot(hfig{plotnum+2},1,1:numtrials,'color',[0 0 0]);
set(hfig{plotnum+2},'box','off')
set(hfig{plotnum+2},'xtick',[],'xcolor',[1 1 1]); pause(0.1);
set(hfig{plotnum+2},'ytick',[],'ycolor',[1 1 1]);  pause(0.1);
xvert=xlabel(hfig{plotnum+2},'big reward side','color',[0 0 0]); 
set(xvert,'units','pixels');
xvertpos=get(xvert,'position');

%xvertpos=getpixelposition(xvert);
set(xvert,'Units','Pixel','Position', [xvertpos(1)+5 xvertpos(2)+colorsize(2)+50],'rotation',90);
hold(hfig{plotnum+2},'on')
for ii=1:numtrials
    if vdata2b(ii)==0
        text(hfig{plotnum+2},1,ii,'\leftarrow','fontsize',8)
    else
        text(hfig{plotnum+2},1,ii,'\rightarrow','fontsize',8)
    end
end
ylim(hfig{plotnum+2},[1 numtrials+1]);
set(hfig{plotnum+2},'YDir','reverse') 

hold(hfig{plotnum+2},'off')

%Previous success/failure
vdata3=abs(data.prevsuccess(seltrials)-1); %invert array to get failures
plot(hfig{plotnum+3},1,1:numtrials,'color',[0 0 0]);
set(hfig{plotnum+3},'box','off')
set(hfig{plotnum+3},'xtick',[],'xcolor',[1 1 1]); pause(0.1);
set(hfig{plotnum+3},'ytick',[],'ycolor',[1 1 1]);  pause(0.1);
xvert2=xlabel(hfig{plotnum+3},'prev fail','color',[0 0 0]); 
set(xvert2,'units','pixels');
xvertpos2=get(xvert2,'position');

%xvertpos2=getpixelposition(xvert2);
set(xvert2,'Units','Pixel','Position', [xvertpos2(1)+5 xvertpos2(2)+colorsize(2)+50],'rotation',90);
hold(hfig{plotnum+3},'on')
for ii=1:numtrials
    if vdata3(ii)==1
        text(hfig{plotnum+3},1,ii,'x','fontsize',8)
    end
end
ylim(hfig{plotnum+3},[1 numtrials+1]);
set(hfig{plotnum+3},'YDir','reverse') 

hold(hfig{plotnum+3},'off')

%Previous rewarded side
vdata4=data.prevside(seltrials);
plot(hfig{plotnum+4},1,1:numtrials,'color',[0 0 0]);
set(hfig{plotnum+4},'box','off')
set(hfig{plotnum+4},'xtick',[],'xcolor',[1 1 1]); pause(0.1);
set(hfig{plotnum+4},'ytick',[],'ycolor',[1 1 1]);  pause(0.1);
xvert3=xlabel(hfig{plotnum+4},'prev side','color',[0 0 0]); 
set(xvert3,'units','pixels');
xvertpos3=get(xvert3,'position');

%xvertpos3=getpixelposition(xvert3);
set(xvert3,'Units','Pixel','Position', [xvertpos3(1)+5 xvertpos3(2)+colorsize(2)+50],'rotation',90);
hold(hfig{plotnum+4},'on')
switchtrials=[];
for ii=1:numtrials
    markc=[0 0 0];
    if vdata4(ii)~=vdata2(ii)
        %if different direction than current trial change color
        markc=[1 0 0];
        if ~isempty(switchtrials)
            switchtrials=[switchtrials ii];
        else 
            switchtrials=ii;
        end
    end
    if vdata4(ii)==0
        text(hfig{plotnum+4},1,ii,'\leftarrow','fontsize',8,'color',markc)
    else
        text(hfig{plotnum+4},1,ii,'\rightarrow','fontsize',8,'color',markc)
    end
end
ylim(hfig{plotnum+4},[1 numtrials+1]);
set(hfig{plotnum+4},'YDir','reverse') 

hold(hfig{plotnum+4},'off')

%Next rewarded side
vdata5=data.nextside(seltrials);
plot(hfig{plotnum+5},1,1:numtrials,'color',[0 0 0]);
set(hfig{plotnum+5},'box','off')
set(hfig{plotnum+5},'xtick',[],'xcolor',[1 1 1]); pause(0.1);
set(hfig{plotnum+5},'ytick',[],'ycolor',[1 1 1]);  pause(0.1);
xvert3=xlabel(hfig{plotnum+5},'next side','color',[0 0 0]); 
set(xvert3,'units','pixels');
xvertpos3=get(xvert3,'position');

%xvertpos3=getpixelposition(xvert3);
set(xvert3,'Units','Pixel','Position', [xvertpos3(1)+5 xvertpos3(2)+colorsize(2)+50],'rotation',90);
hold(hfig{plotnum+5},'on')
switchnexttrials=[];
for ii=1:numtrials
    markc=[0 0 0];
    if vdata5(ii)~=vdata2(ii)
        %if different direction than current trial change color
        markc=[1 0 0];
        if ~isempty(switchnexttrials)
            switchnexttrials=[switchnexttrials ii];
        else 
            switchnexttrials=ii;
        end
    end
    if vdata5(ii)==0
        text(hfig{plotnum+5},1,ii,'\leftarrow','fontsize',8,'color',markc)
    else
        text(hfig{plotnum+5},1,ii,'\rightarrow','fontsize',8,'color',markc)
    end
end
ylim(hfig{plotnum+5},[1 numtrials+1]);
set(hfig{plotnum+5},'YDir','reverse')        %flip y-axis values so first trial on top
hold(hfig{plotnum+5},'off')

%Next success/failure
vdata6=abs(data.nextsuccess(seltrials)-1); %array to get failures
plot(hfig{plotnum+6},1,1:numtrials,'color',[0 0 0]);
set(hfig{plotnum+6},'box','off')
set(hfig{plotnum+6},'xtick',[],'xcolor',[1 1 1]); pause(0.1);
set(hfig{plotnum+6},'ytick',[],'ycolor',[1 1 1]);  pause(0.1);
xvert5=xlabel(hfig{plotnum+6},'next fail','color',[0 0 0]); 
set(xvert5,'units','pixels');
xvertpos2=get(xvert5,'position');

%xvertpos2=getpixelposition(xvert5);
set(xvert5,'Units','Pixel','Position', [xvertpos2(1)+5 xvertpos2(2)+colorsize(2)+50],'rotation',90);
hold(hfig{plotnum+6},'on')
for ii=1:numtrials
    if vdata6(ii)==1
        text(hfig{plotnum+6},1,ii,'x','fontsize',8)
    end
end
ylim(hfig{plotnum+6},[1 numtrials+1]);
set(hfig{plotnum+6},'YDir','reverse')        %flip y-axis values so first trial on top
hold(hfig{plotnum+6},'off')

pos1=getpixelposition(hfig{1});
%%%draw rectangles to visualize individual trial separations
set(hfig{plotnum+7}, 'Units','Pixels','Position',  [origpos(1) pos1(2) vertplotwidth(1)+vertplotwidth2*6+10 colorsize(2)]);
%x=[0 1]; x=repmat(x,1,numtrials); y=sort(repmat(1:numtrials,1,2));
x=[0 0;  1 1]; x1=repmat(x,[floor(numtrials/2),1]);
cla(hfig{plotnum+7});
delete(findobj(hfig{plotnum+7},'type','line'));         %delete previous markers trial #'s
p1=image(x1,'parent',hfig{plotnum+7},'cdatamapping','scaled');                 %make tiled rows every other trial
alpha(p1,0.2);
ylim(hfig{plotnum+7},[0.5 numtrials+.5])
cdata=repmat([1 1 1;0.7 0.7 0.7],ceil(numtrials/2),1);
   % 

colormap(hfig{plotnum+7},cdata);       %black/white color map tiled rows

hold(hfig{plotnum+7},'on')
set(hfig{plotnum+7},'color','none');        %transparent background
set(hfig{plotnum+7},'box','off');
set(hfig{plotnum+7},'xtick',[],'xcolor',[1 1 1]); pause(0.1);
set(hfig{plotnum+7},'ytick',[],'ycolor',[1 1 1]);  pause(0.1);
%{
%yyaxis right;      %not available in 2013
[aa,p1,p2]=plotyy(1,.5:1:numtrials+.5,1,.5:1:numtrials+.5,'parent',hfig{plotnum+7});

set(p2,'marker','none');
ylim(aa(2),[0.5 numtrials+.5]);
yticks2=1:numtrials;
yticks2label=seltrials;
set(aa(2),'ytick',yticks2,'ycolor',[0 0 0]); 
set(aa(2),'tickdir','out');
set(aa(2),'yticklabel',yticks2label)
set(aa(2),'fontsize',7)
hold(aa(2),'off')
set(aa(2),'YDir','reverse')        %flip y-axis values so first trial on top

%}
%if reward side is left (0) get reaction times for both sides separate
rts.left=vdata(vdata2==0);
rts.right=vdata(vdata2==1);

end