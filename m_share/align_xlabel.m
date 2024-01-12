function align_xlabel(hxlabel,posbot)

% function align_xlabel(hxlabel,posbot)
% posbot is optional (in normalized units, bottommost position for the
% Xlabel)


for i=1:length(hxlabel)
set(hxlabel(i),'units','normalized')
hpos(i,:)=get(hxlabel(i),'position');
end
miny=min(hpos(:,2));
if nargin>1
    miny=posbot;
end

for i=1:length(hxlabel)
hposnew(i,:)=hpos(i,:); %hpos1new(1)=-0.07 %here you put what you want
hposnew(i,2)=miny;
set(hxlabel(i),'position',hposnew(i,:),'verticalalignment','bottom');
end


