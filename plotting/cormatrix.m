function cor=cormatrix(d1,d2,types,varargin)
%2/19/2019, save d1site & d2site as char, not as cell (prev d1sites were
%all cells), d2 sites were ok
numchs=size(d1,2);
count=0;
while 1
    count=count+1;
    if ~isempty(d1{count})
        break
    end
end
firstch=count;
argnum=1;
sametype=0;
sitelab=0;
logtype=0;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'sametype'
            sametype=1;
        case 'sitelabels'
            sitelab=1;
        case 'log'
            logtype=1;
    end
    argnum=argnum+1;
end
allfields=fieldnames(d1{firstch});

rdata={};
pdata={};
corrdata={};
cor={};
for ich=firstch:numchs
            corcount=1;
for ifield=1:length(types)
    if ~ismember(types{ifield},allfields)
        error('no such field type');
    end
    fieldtype=types{ifield};
    if contains(fieldtype,'diff')
        continue
        %diff types are same for da signals, but not for lfp (iifields)
    end
    %for each channel in d1 correlate to all channels in d2
    if ~isempty(d1{ich})
        d1data=getfield(d1{ich},types{ifield});
        if size(d1data,2)>size(d1data,1)
                d1data=d1data';
        end
        d1trials=d1{ich}.trialnums;
        %get correlation of da ch against serial sequence of trials
        rs=nan;
        ps=nan;
        if ~isempty(d1data)
        [rs,ps]=corr(d1data,d1{ich}.trialnums','rows','complete');
        end
        if ps<=0.05
        %if correlated wiht progression of task, then need to be
        %careful with interpretation
        cor{ich}.corrdata(corcount).type1=['task progression-' types{ifield}];
        cor{ich}.corrdata(corcount).type2='na';
        cor{ich}.corrdata(corcount).d1ch=ich;
        d1site=getfield(d1{ich},'site');
        if iscell(d1site)
            d1site=d1site{:};
        end
        cor{ich}.corrdata(corcount).d1site=d1site;
        cor{ich}.corrdata(corcount).d2ch=0;
        cor{ich}.corrdata(corcount).d2site='na';
        cor{ich}.corrdata(corcount).r=rs;
        cor{ich}.corrdata(corcount).p=ps;
        corcount=corcount+1;
        end
        for iich=1:size(d2,2)
                if isempty(d2{iich})
                    continue
                end
            %get data for each channel in d2 data
            if sametype==0
                scantypes=1:length(types);
                %look at all other types for cross-correlations
            else
                %scantypes=ifield;
                %look at only same type field
                %or with same name wihin ie rewimpeak rewwin etc.
                %by looking at first 3 letters of type
                scantypes=find(contains(types,types{ifield}(1:3))==1);
            end
            for iifield=scantypes

                d2data=getfield(d2{iich},types{iifield});
                d2trials=d2{iich}.trialnums;
                if size(d2data,2)>size(d2data,1)
                    d2data=d2data';
                end
                if logtype==1
                    d2data=log10(d2data);
                end
                if ich==firstch && ifield==1
                %get correlation lfp ch against serial sequence of trials
                rs=nan;
                ps=nan;
                if ~isempty(d2data)
                [rs,ps]=corr(d2data,d2{iich}.trialnums','rows','complete');
                end
                if ps<=0.05
                %if correlated wiht progression of task, then need to be
                %careful with interpretation
                cor{ich}.corrdata(corcount).type1='na';
                cor{ich}.corrdata(corcount).type2=['task progression-' types{iifield}];
                cor{ich}.corrdata(corcount).d1ch=0;
                cor{ich}.corrdata(corcount).d1site='na';
                cor{ich}.corrdata(corcount).d2ch=iich;
                d2site=getfield(d2{iich},'site');
                if iscell(d2site)
                    d2site=d2site{:};
                end
                cor{ich}.corrdata(corcount).d2site=d2site;
                cor{ich}.corrdata(corcount).r=rs;
                cor{ich}.corrdata(corcount).p=ps;
                corcount=corcount+1;
                end
                end
                %use trialnums to find same trial data since added new
                %ch-specific good trials 05/2019
                commontrials=intersect(d1trials,d2trials);
                d1ids=find(ismember(d1trials,commontrials));
                d2ids=find(ismember(d2trials,commontrials));
                rs=nan;
                ps=nan;
                if ~isempty(d1data(d1ids)) && ~isempty(d2data(d2ids))
                    [rs,ps]=corr(d1data(d1ids),d2data(d2ids),'rows','complete');
                end

                if ps<=0.05
                    cor{ich}.corrdata(corcount).type1=types{ifield};
                    cor{ich}.corrdata(corcount).type2=types{iifield};
                    cor{ich}.corrdata(corcount).d1ch=ich;
                    d1site=getfield(d1{ich},'site');
                    if iscell(d1site)
                        d1site=d1site{:};
                    end
                    cor{ich}.corrdata(corcount).d1site=d1site;
                    cor{ich}.corrdata(corcount).d2ch=iich;
                    d2site=getfield(d2{iich},'site');
                    if iscell(d2site)
                        d2site=d2site{:};
                    end
                    cor{ich}.corrdata(corcount).d2site=d2site;
                    cor{ich}.corrdata(corcount).r=rs;
                    cor{ich}.corrdata(corcount).p=ps;
                    corcount=corcount+1;
                end
                cor{ich}.rdata=setfield(rdata,{ich},types{ifield},...
                    types{iifield},{iich},rs);
                cor{ich}.pdata=setfield(pdata,{ich},types{ifield},...
                    types{iifield},{iich},ps);
            end

        end
    else
        cor{ich}={};
    end
end
    %put 'task progression' types at top
    %sort correlation matrix by ascending p
if ~isempty(cor{ich}) && isfield(cor{ich},'corrdata')
ps=[cor{ich}.corrdata.p];
[a,b]=sort(ps);
corsort=cor{ich}.corrdata(b);
cor{ich}.corrdata=corsort;
aa={cor{ich}.corrdata(:).type1};
ab=contains(aa,'task');
bb=find(ab==1);
aa={cor{ich}.corrdata(:).type2};
abc=contains(aa,'task');
bbc=find(abc==1);
bb=unique([bb bbc]);
ab=ab+abc;
amat=cor{ich}.corrdata(bb);
cc=find(ab==0);
amat2=cor{ich}.corrdata(cc);
amat3=[amat amat2];
cor{ich}.corrdata=amat3;
end
end
  
%{
%sort correlation matrix by ascending p
ps=[corrdata.p];
[a,b]=sort(ps);
corsort=corrdata(b);
corrdata=corsort;
%}
end
        
        
        