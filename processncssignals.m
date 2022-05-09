function processed=processncssignals(samples,samplerate,idgroups,filtgroups,varargin)
%filter & envelope signals
%samples is ch in x and samples in y
%idcats is struct array where each cell has channels for each category
%ie. idcats{1}=(1:10), means first 10 chs correspond filtcats{1}=[17 34];
argnum=1;
envgroups=[];
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'envgroups'
            argnum=argnum+1;
            envgroups=varargin{argnum};
    end
    argnum=argnum+1;
end
if size(samples,1)>size(samples,2)
    %want ch to be in row
    samples=samples';
end
        
ichidx=1;
for id=1:length(idgroups)
    samplesgroup=samples(idgroups{id},:);
    filt=filtgroups{id};
    for ich=1:size(samplesgroup,1)
        if ~isempty(filt)
            %if filt band indicated
            processed(ichidx,:)=filterLFP(samplesgroup(ich,:),samplerate,filt);
        else
            %retain original input signal
            processed(ichidx,:)=samplesgroup(ich,:);
        end
        if ~isempty(envgroups)
            if envgroups(id)==1
                envwin=round(samplerate*.5/mean(filt)); %envelop over 1/2 cycle of targetd osc 
                squaredsig=processed(ichidx,:).^2;   %get power V^2
                processed(ichidx,:)=smoothwin(squaredsig,envwin);   %smoothing
            end
        end
        ichidx=ichidx+1;
    end
end
    
end