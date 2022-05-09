function togglee(source,event)
%erase detected signals if overlap with badids (button m and manually clik)
global plotParam processed parameters hgui
pressed = source.Value;     %checked?  
selch=plotParam.selch;
ovwin=parameters.samplerate*2;      %2 s overlap remove
if pressed==1 && ~isempty(processed.badids)
    badids=processed.badids;
    for ii=1:length(selch)
        if ii<=length(processed.detected)
            %if detected signals exist for ch
            if ~isempty(processed.detected{selch(ii)})
                da=processed.detected{selch(ii)};
                %only if detected signals exist for this channel
                tids=processed.detected{selch(ii)}.tids;
                tracelength=size(tids,2);
                pid=median(1:tracelength);
                dids=tids(:,pid);       %get id just at peak ts's for all det
            for id=1:length(badids)
                %scroll all bad ids identified and find overlapping time
                %periods of detected peaks within 1 s remove
            ovids=intersect(dids,badids(id)-ovwin:badids(id)+ovwin);
                if ~isempty(ovids)
                    %if  identified, remove
                    %RAther than remove, just identify as bad id & update
                    %selected ids in xcov variable
                    %make da trace nan
                    nib=find(ismember(dids,ovids)==1); 
                    %da.maxTS(nib)=[];
                    da.datrace(nib,:)=nan(1,size(da.datrace,2));   
                    for ilfp=1:length(da.lfptrace)
                        %da.lfptrace{ilfp}(nib,:)=[];
                        %da.nlxtrace{ilfp}(nib,:)=[];
                        xib=find(ismember(da.xcov{ilfp}.sel,nib)==1); 
                        if ~isempty(xib)
                            %remove as selected trial, if was selected
                            da.xcov{ilfp}.sel(xib)=[];
                        end
                    end
                    if id==1
                        da.badpeaks=nib;
                    else
                        da.badpeaks=[da.badpeaks nib];
                    end
                end
            end
            da.badids=badids;
            processed.detected{selch(ii)}=da;   %replace original processed variable
            end
        end
    end
                   
else

    disp('no bad ids to remove detected peaks with');
end
set(hgui.erase,'Value',0);
end
