color=hsv;
color=downsample(color,round(length(color)/(numch)));
if numch>4
    lightcolors1=color(4,:);
    lightcolors2=color(5,:);
    color(4,:)=[lightcolors1(1) lightcolors1(2)-.2 lightcolors1(3)];
    color(5,:)=[lightcolors2(1)-.1 lightcolors2(2)-.2 lightcolors2(3)];
    color1=[.3 0 0]; color2=[.6 0 0];
    colorCV=[color1; color2; color];
else
colorCV=color;
end