function data=gettcov(sessnums,xinfos,datm,plotparam,twin,ttypes,varargin)
twin=[-3 4];
outdataall=[];
savepath=plotparam.savepath;
savepath=[savepath filesep 'tcov' filesep];
if ~isdir(savepath)
    mkdir(savepath);
end
argnum=1;
lfpsites={};
dasites={};
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'lfps'
            %specify lfpsite
            argnum=argnum+1;
            lfpsites=varargin{argnum};
        case 'dasites'
            argnum=argnum+1;
            dasites=varargin{argnum};
    end
    argnum=argnum+1;
end

for ise=1:length(sessnums)
cursess=sessnums(ise)
if isempty(lfpsites)
outdatatemp=plotmultxtraces(cursess,xinfos,plotparam,'win',twin,'event','targ','ttypes',ttypes,'datm',datm,'sort','targpeak','tcov');
else
    outdatatemp=plotmultxtraces(cursess,xinfos,plotparam,'win',twin,'event','targ','ttypes',ttypes,'datm',datm,'sort','targpeak','tcov','lfps',lfpsites,'dasites',dasites);
end
if ise==1
outdataall=outdatatemp;
else
outdataall=[outdataall outdatatemp];
end
close all;
end

sumx=[];
stdv=[];
pmask=[];
for ix=1:length(outdataall)
if ix==1
sumx=outdataall(1).xvar;
else
sumx=outdataall(ix).xvar+sumx;
end
end
conditions=[];
for it=1:length(ttypes)
    if it==1
    conditions=ttypes{1};
    else
        conditions=[conditions '_' ttypes{it}];
    end
end
%sig mask
for ix=1:size(outdataall(1).xvar,1)
    %rows
    for iy=1:size(outdataall(1).xvar,2)
        %columns
        pvals=[];
        for xx=1:length(outdataall)
            %if sum(sum(outdataall(xx).p<.05))>.03*71*71
                pvals(xx)=outdataall(xx).p(ix,iy);
            %else
            %    disp([num2str(outdataall(xx).sessnum) 'not sig']);
            %    pvals(xx)=nan;
            %end
        end
                pmask(ix,iy)=sum(pvals);
                pmask(ix,iy)=dg_chi2test2([sum(pvals>=0.05) sum(pvals>0)]);%count test...
            pmask(ix,iy)=-2*sum(log(pvals)); %fisher's test...
    end
end

figpos=[50,50,1000,800];
axsize=[600 600];
figsess=figure('visible','off');     %figure for each channel
if ispc
figsess=figure('visible','on');     %figure for each channel
end
set(figsess,'position',figpos,'color',[1 1 1]);
set(0,'CurrentFigure',figsess);    %set figure handle to current figure
axa=axes;

avgx=sumx./length(outdataall);
imagedata=imagesc(axa,avgx)
 artTime=(pmask==1);   %find artifact points (nan periods)
        maskGray=artTime*.01;            %make gray rather than white default by making non-zero
        artTime=artTime+maskGray;
        set(imagedata, 'AlphaData', artTime); 
        artdata=maskGray+avgx;
        im2=imagesc(axa,artdata);
interval=1;
rate=10;
xticklabels=min(twin):interval:max(twin);
xticklabels=round(xticklabels.*rate)./rate;
xticklabels=num2str(xticklabels');
xticks=1:round(interval*rate):size(avgx,2);
set(axa,'xtick',xticks,'xticklabel',xticklabels);
set(axa,'xlim',[1 size(avgx,2)]);
set(axa,'tickdir','out','box','off')
xlabel(axa,'DA time relative to target cue (s)')
set(axa,'ytick',xticks,'yticklabel',xticklabels);
set(axa,'ylim',[1 size(avgx,2)]);
ylabel(axa,'\beta time relative to target cue (s)','interpreter','tex')
set(axa,'ylim',[1 size(avgx,1)]);
origpos=getpixelposition(axa);      %get original position        
h1=colorbar(axa,'eastoutside');  
set(h1,'units','pixels');
set(axa,'Units','Pixels','Position', [origpos(1) origpos(2) axsize(1) axsize(2)]);
%clabel('average r')
%set(h1,'position',[axpos{ip}(1) 30 100 10]);
set(findall(axa,'-property','FontSize'),'FontSize',20)
aa=title(h1,'correlation coefficient','interpreter','tex')
%set([-0.04 0.04])
set(aa,'rotation',270)
set(aa,'units','pixels')
cpos=get(aa,'position');
set(aa,'position',[cpos(1)+100 cpos(2)/2]);
ht=title(axa,'Session averaged covariance of dopamine as a function of time and beta as a function of time','fontsize',10)
set(ht,'units','pixels')
hpos=get(ht,'position');
set(ht,'position',[hpos(1) cpos(2)+20]);
set(axa,'xlim',[31 size(avgx,2)]);
set(axa,'ylim',[31 size(avgx,2)]);
set(axa,'YDir','normal');    %flip up down
xmap=brewermap(100,'RdBu');
colormap(xmap)
caxis([-0.04 0.04])

savename=[savepath  'timerelcovariance_allsessions' '_' 'targalign' '_' 'paired beta-da sites_' conditions];
if length(sessnums)<=1
    savename=[savename '_' num2str(sessnums(1))];
end
savefig(figsess,savename);
saveas(figsess,savename,'jpg')
print(figsess,savename,'-painters','-depsc');
data=outdataall;
