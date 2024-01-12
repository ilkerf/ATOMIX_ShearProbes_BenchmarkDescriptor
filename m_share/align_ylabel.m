function align_ylabel(hylabel,posleft)

% function align_ylabel(hylabel,posleft)
% posleft is optional (in normalized units, left-most position for the
% ylabel)


for i=1:length(hylabel)
set(hylabel(i),'units','normalized')
hpos(i,:)=get(hylabel(i),'position');
end
miny=min(hpos(:,1));
if nargin>1
    miny=posleft;
end

for i=1:length(hylabel)
hposnew(i,:)=hpos(i,:); %hpos1new(1)=-0.07 %here you put what you want
hposnew(i,1)=miny;
set(hylabel(i),'position',hposnew(i,:));
end


