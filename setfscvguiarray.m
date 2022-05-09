function setfscvguiarray(fh,numcolor)
%set up figure & individual plots
global plotParam parameters hgui
clf(fh);
set(0,'CurrentFigure',fh);    %set figure handle to current figure
hfig={};
hpos={};
figpos=[50,50,1500,800];
hgui.loadedcsc=plotParam.cscNames;
plotsize=plotParam.colorplotsize;
widen=plotParam.widen;
%figpos=get(fh,'position');
set(fh,'position',figpos)

hgui.hf=fh;
hgui.txtpad=25;
%set up plotting axes
for ifig=1:numcolor
    %subplots 1, 4, 7, 10 are color plots
    if ifig<=4
    hgui.hfig{ifig}=subplot(5,5,(ifig-1)*5+1);     %1, 6, 11, 16
    end
    if ifig<=8 && ifig>4
    hgui.hfig{ifig}=subplot(5,5,(ifig-4-1)*5+2);     %2, 7, 12, 17      
    end
    if ifig<=12 && ifig>8
    hgui.hfig{ifig}=subplot(5,5,(ifig-8-1)*5+3);     %3, 8, 13, 18        
    end
    if ifig<=16 && ifig>12
    hgui.hfig{ifig}=subplot(5,5,(ifig-12-1)*5+4);     %4, 9, 14, 19          
    end
end
hgui.titletext=subplot(5,5,5);
hgui.itplot=subplot(5,5,21);
hgui.itplotx{1}=subplot(5,5,22);
hgui.itplotx{2}=subplot(5,5,23);
hgui.itplotx{3}=subplot(5,5,24);
hgui.cv=subplot(5,5,10);         
hgui.closeup{1}={};         
hgui.data=subplot(5,5,15);         
hgui.closeup{2}={};        
hgui.closeup{3}={};    
hgui.fftplot={};       

for ifig=1:numcolor
    hpos{ifig}=getpixelposition(hgui.hfig{ifig});
end
hpostitle=getpixelposition(hgui.titletext);
hposit=getpixelposition(hgui.itplot);
hposcv=getpixelposition(hgui.cv);
hposdata=getpixelposition(hgui.data);

hposcloseup{1}=getpixelposition(hgui.closeup{1});
hposcloseup{2}=getpixelposition(hgui.closeup{2});
hposcloseup{3}=getpixelposition(hgui.closeup{3});

%resize color plots stretch
margins=25;
titlewidth=450; titleheight=10;
xoff=100;
xmarg=125;
yoff=8;
yinitialoff=75;
countb=0;
ymarg=40;
inity=hpos{1}(2);

%color plots hgui.hfig{1:4}
for ifig=1:numcolor
    set(hgui.hfig{ifig}, 'Units','Pixels');
    set(hgui.hfig{ifig},'Position', ...
        [xoff+(ceil(ifig/4)-1)*plotsize(1)+margins*(ceil(ifig/4)-1) figpos(4)-yinitialoff-plotsize(2)*(ifig-(ceil(ifig/4)-1)*4-1)-yinitialoff-margins*(ifig-(ceil(ifig/4)-1)*4-1) ...
        plotsize(1) plotsize(2)]);
    hpos{ifig}=getpixelposition(hgui.hfig{ifig});
    %{
    if ~ismember(ifig,plotParam.selch)
        set(hgui.hfig{ifig},'ytick',[],'xtick',[],'visible','off')
    end
    %}
end

%pca vs time plot hgui.itplot
set(hgui.itplot, 'Units','Pixels','Position',  ...
    [xoff hpos{4}(2)-yoff*(4)-yinitialoff*1.75 ...
    plotsize(1) plotsize(2)]);
itpos=get(hgui.itplot,'position');
set(hgui.itplotx{1}, 'Units','Pixels','Position',  ...
    [hpos{8}(1) hpos{8}(2)-yoff*(4)-yinitialoff*1.75 ...
    plotsize(1) plotsize(2)]);
itposx1=get(hgui.itplotx{1},'position');
set(hgui.itplotx{2}, 'Units','Pixels','Position',  ...
    [hpos{12}(1) hpos{12}(2)-yoff*(4)-yinitialoff*1.75 ...
    plotsize(1) plotsize(2)]);
itposx2=get(hgui.itplotx{2},'position');
set(hgui.itplotx{3}, 'Units','Pixels','Position',  ...
    [hpos{16}(1) hpos{16}(2)-yoff*(4)-yinitialoff*1.75 ...
    plotsize(1) plotsize(2)]);
itposx3=get(hgui.itplotx{3},'position');

%cv plot hgui.cv
set(hgui.cv, 'Units','Pixels','Position',  ...
    [hpos{16}(1)+plotsize(1)+xoff-60 figpos(4)-250 plotsize(1) plotsize(2)]);
pos6=getpixelposition(hgui.cv);

%text data plot hgui.data
set(hgui.data, 'Units','Pixels','Position',  ...
    [pos6(1) pos6(2)-200 plotsize(1) plotsize(2)]);

%title text
set(hgui.titletext,'Units','Pixels','Position',...
    [figpos(3)-titlewidth figpos(4)-titleheight*2 titlewidth titleheight]);
set(hgui.titletext,'ytick',[],'xtick',[],'visible','off')

%set up buttons
hgui.menu = uicontrol('Style', 'popup',...
   'String', {'menu','load file','load settings',...
   'save settings','select nlx chs',...
   'fft scale',...
   'detect da unbias (current window)','detect da unbias (entire file)',...
   'get xcov da w/ lfp (current window)','get xcov da w/ lfp (entire file)','merge select files',...
    'get da for directory','get da and xcov for directory','load detected da auto dir',...
   'save detected da','load bad ids','save bad ids','save trials to xcel'},...
   'Position', [figpos(3)-100 figpos(4)-80 100 25],...
   'Callback', @menuselect);  
%check plots if in plotParam.selch
checkoff=-30;
for ich=1:numcolor
    plotvalue=0;
    if ismember(ich,plotParam.selch)
        plotvalue=1;
    end
    if ich>4
        checkoff=0;
    end
    hgui.check{ich} = uicontrol('Style', 'checkbox','value',plotvalue,... 
        'string',[num2str(ich)],'Position', [hpos{ich}(1)-20+checkoff ...
        hpos{ich}(2)+plotsize(2)-15 15 15],...
        'backgroundcolor',plotParam.colorFSCV((ich-4*(ceil(ich/4)-1)),:),...
        'callback',@checkbox);      
     checktxt = uicontrol('Style', 'text','string',num2str(ich),...
        'Position', [hpos{ich}(1)-20+checkoff hpos{ich}(2)+plotsize(2)-35 20 15],...
        'backgroundcolor',[1 1 1],...
        'foregroundcolor',plotParam.colorFSCV((ich-4*(ceil(ich/4)-1)),:));  
end

%text inputs (time periods to plot & color max scale)
textt = uicontrol('Style', 'text',...
   'String', 'start / end (s)',...
   'Position', [figpos(3)-185 figpos(4)-130 70 25],...
   'BackgroundColor', [ 1 1 1]);  
hgui.tstart = uicontrol('Style', 'edit',...
   'String', num2str(round(plotParam.t_start/parameters.samplerate)),...
   'Position', [figpos(3)-100 figpos(4)-120 40 25]);  
hgui.tend = uicontrol('Style', 'edit',...
   'String', num2str(round(plotParam.t_end/parameters.samplerate)),...
   'Position', [figpos(3)-50 figpos(4)-120 40 25]);  

hgui.refresh=uicontrol('style','pushbutton','string','refresh',...
    'position',[figpos(3)-100 figpos(4)-200 50 25],'callback',@refreshbutton);

textcm = uicontrol('Style', 'text',...
   'String', 'I max',...
   'Position', [figpos(3)-150 figpos(4)-170 40 25],...
   'BackgroundColor', [ 1 1 1]);  
hgui.cmax = uicontrol('Style', 'edit',...
   'String', num2str(plotParam.cmax),...
   'Position', [figpos(3)-100 figpos(4)-160 25 25]);  

hposit=getpixelposition(hgui.itplot);

%plot bad ids
 hgui.badplot= uicontrol('Style', 'togglebutton','value',0,... 
        'string','m','Position', [hposit(1)-70 ...
        hposit(2)+plotsize(2)-25 20 20],...        
        'callback',@togglem);  
%erase detected peaks based on bad ids identified
 hgui.erase= uicontrol('Style', 'togglebutton','value',0,... 
        'string','E','Position', [hposit(1)-90 ...
        hposit(2)+plotsize(2)-25 20 20],...        
        'callback',@togglee);  
%scroll time windows
 hgui.lscroll= uicontrol('Style', 'pushbutton',... 
        'string','<<','Position', [hposit(1)-90 ...
        hposit(2) 20 20],...        
        'callback',@scroll); 
  hgui.rscroll= uicontrol('Style', 'pushbutton',... 
        'string','>>','Position', [hposit(1)-70 ...
        hposit(2) 20 20],...        
        'callback',@scroll); 
       
%text legend for lfp traces plotted
%updatelegends(hgui);
updatelegends;

%{
cscnames=plotParam.cscNames;
lfpid=plotParam.lfpid;
physid=plotParam.physid;
colortab=plotParam.colormaptab;
txtpad=30;
for ii=1:length(lfpid)
hgui.labelslfp{ii}= uicontrol('Style', 'text',...
   'String', cscnames{lfpid(ii)},'fontangle','italic',...
   'Position', [figpos(3)-60 ...
   hposcloseup{1}(2)-50-txtpad*(ii-1) 60 25],...
   'BackgroundColor', [ 1 1 1],'foregroundcolor',colortab(ii,:));
end
for ii=1:length(physid)
hgui.labelslfp{ii}= uicontrol('Style', 'text',...
   'String', cscnames{physid(ii)},'fontangle','italic',...
   'Position', [figpos(3)-60 ...
   hposcloseup{1}(2)-450-txtpad*(ii-1) 60 25],...
   'BackgroundColor', [ 1 1 1],'foregroundcolor',colortab(ii,:));
end
hgui.txtpad=txtpad;
%}
colormap(plotParam.map)

hgui.itplotpos{1}=get(hgui.itplot,'position');
hgui.itplotpos{2}=get(hgui.itplotx{1},'position');
hgui.itplotpos{3}=get(hgui.itplotx{2},'position');
hgui.itplotpos{4}=get(hgui.itplotx{3},'position');


end
