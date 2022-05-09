function grpdata=plothistogroups(response,rvar,varargin)
argnum=1;
plothist=0;
fontsize=12;
units=[];
varnames=[];
savename=[];
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'plot'
            plothist=1;
        case 'responsename'
            argnum=argnum+1;
            units=varargin{argnum};
        case 'response'
            argnum=argnum+1;
            %get response variable from data provided
            %ie. getfield(datmall{4},'targabspeak')
            targfield=varargin{argnum};
            units=targfield;
            response=getfield(response,targfield);
        case 'varnames'
            argnum=argnum+1;
            varnames=varargin{argnum};
        case 'savename'
            argnum=argnum+1;
            savename=varargin{argnum};
    end
    argnum=argnum+1;
end
if isempty(units)
    units='response';
end
if isempty(varnames)
    for ix=1:length(rvar)
        varnames{ix}='';
    end            
end
grpdata=[];
binwidth=2.5;
binmin=min(response);
binmax=max(response);
binmin=-binmax;
edges=binmin:binwidth:binmax;
itypes1=unique(rvar{1});
itypes2=unique(rvar{2});



count=1;
for it=1:length(itypes1)    
    for iit=1:length(itypes2)
        grptids=find(strcmp(rvar{1},itypes1{it}) & strcmp(rvar{2},itypes2{iit}));
        grpdata.response{count}=response(grptids);
        grpdata.ratioposneg{count}=[{grptids(response(grptids)>0)} {grptids(response(grptids)<0)}];
        grpdata.variable{count}={itypes1{it} itypes2{iit}};
        count=count+1;
    end
end

if plothist
    fighist=figure('position',[50 50 1500 800],'color',[1 1 1]);
    ax={};
    iax=1;
    for it=1:length(itypes1)
        for iit=1:length(itypes2)
            ax{it,iit}=subplot(length(itypes1),length(itypes2),iax);
            histogram(ax{it,iit},grpdata.response{iax},edges);
            set(ax{it,iit},'units','pixels','box','off');
            titletext={itypes1{it} itypes2{iit}};
            sizax=get(ax{it,iit},'position');
            text(ax{it,iit},20,sizax(4),{[varnames{1} ' ' itypes1{it}] [varnames{2} ' ' itypes2{iit}]},'units','pixels','fontsize',fontsize);
            if it==length(itypes1)
                xlabel(ax{it,iit},units);
            end
            xlim(ax{it,iit},[binmin binmax]);
            ylabel(ax{it,iit},'counts');
            iax=iax+1;
        end
    end
    if ~isempty(savename)
        savefig(fighist,savename);
        saveas(fighist,savename,'tif')
    end
end

end