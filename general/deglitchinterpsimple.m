function deglitched=deglitchinterpsimple(data,thres,width)
deglitched=data;
for ix=1:size(data,1)
    diffdata=diff(data(ix,:));
    glitchpt=find(abs(diffdata)>abs(nanstd(diffdata))*thres)+1;
    if length(glitchpt)>1
        newglitchpts=glitchpt(1);
        count=1;
        for ip=2:length(glitchpt)
            if glitchpt(ip-1)~=glitchpt(ip)-1
                count=count+1;
                newglitchpts(count)=glitchpt(ip);
            end
        end
        glitchpt=newglitchpts;
    end
    for ip=1:length(glitchpt)
    xqs=glitchpt(ip)-width:glitchpt(ip)+width;
    xqs=xqs(xqs>0 & xqs<=length(deglitched(ix,:)));
    deglitched(ix,xqs)=nan;
    xqs2=[xqs(1)-1 xqs xqs(end)+1];
    xqs2=xqs2(xqs2>0 & xqs2<=length(deglitched(ix,:)));
    datainterp=interp1([xqs2(1) xqs2(end)],deglitched(ix,[xqs2(1) xqs2(end)]),xqs2);
    deglitched(ix,xqs2)=datainterp;
    end
end

end
