function [htext,hnewax]=sublabel(varargin)
% function [htext]=sublabel
%  sublabel labels sub-plots from haxes in order as "a", "b",
%  "c".  By default it puts them 12 points from the upper left corner.
%
% can also be called:
%   [htext]=sublabel(ax,dx,dy)
% where ax are the axes handles to be labeled.  dx and dy are offsets from
% the upper left corner (in points).
%
%%ILKER ADDED: loc (location) option
% can also be called:
%   [htext]=sublabel(ax,dx,dy,loc)
% loc = 1  top left (default)
% loc = 2  bottom left
% loc = 3  bottom right
% loc = 4  top right
% Ilker also revised such that loc can be different for each axis
%
% Can also be called:
%   [htext]=sublabel(ax,dx,dy,loc,varargin);
%   [htext]=sublabel(varargin);
%   [htext]=sublabel(ax,varargin);
%   [htext]=sublabel(ax,dx,varargin);
% where varargin is any pair of properties valid for text objects
%  eg:
%    >> haxes(1) = subplot(2,1,1);
%    >> plot(1:10);
%    >> haxes(2) = subplot(2,1,2);
%    >> plot(1:2:20);
%    >> htext=sublabel(haxes);
% Or:
%    >> htext =
%    sublabel(haxes,'fontsize',8,'fontweight','bold','backgroundcolor','g');
%    >> htext =
%    sublabel(haxes,24,36,'fontsize',8,'fontweight','bold','backgroundcolor','g');

% $Revision: 1.1 $ $Date: 2005/04/04 15:36:24 $ $Author: jklymak $
% J. Klymak.  August 8, 2000...

fsize = 10; % default fontsize, points....
fname = 'arial'; % font...
fontw='b';
haxes=flipud(datachildren(gcf));
toffset=fsize+.5;
roffset=fsize-2;%+1;
loc= ones(size(haxes));
num = 0;
while ~isempty(varargin) & ~isstr(varargin{1})
    num = num+1;
    if ~isempty(varargin{1})
        if num==1
            haxes=varargin{1};
        elseif num==2
            toffset=varargin{1}-3;
        elseif num==3
            roffset=varargin{1}-6;
        elseif num==4
            loc=varargin{1};
        else
            error('invalid parameter/value pair.');
        end
    end
    varargin=varargin(2:end);
end
% toffset=fsize;
% roffset=fsize;

if length(loc)<2
    loc=loc.*size(haxes);
end
loc(loc==0)=NaN;
for i=1:length(haxes)
    uni = get(haxes(i),'units');
    hnewax(i)=axes('units',uni,'pos',get(haxes(i),'pos'));
    set(hnewax(i),'visible','off','units','points');
    poss = get(hnewax(i),'pos'); % this is how big the axis is in points...
    axis([poss(1) poss(1)+poss(3) poss(2) poss(2)+poss(4)]);
    x = poss(1)+roffset;
    y = poss(2)+toffset;
    switch loc(i)
        case 1 %top left
            set(hnewax(i),'ydir','rev');
        case 2 %bottom left
            set(hnewax(i),'ydir','n');
            set(hnewax(i),'xdir','n');
        case 3 %bottom right
            set(hnewax(i),'xdir','rev');
            set(hnewax(i),'ydir','n');
       case 4 %top right
            set(hnewax(i),'xdir','rev');
            set(hnewax(i),'ydir','rev');

    end

% Uncomment / Comment for a b c or (a) (b) (c)
    % htext(i)=text(x,y,0,sprintf('%c','(',['a'+i-1],')'),'fontname', ...
    %     fname,'fontsize',fsize,'fontweight',fontw,'backgroundcolor','w',varargin{:});
    htext(i)=text(x,y,0,sprintf('%c',['a'+i-1]),'fontname', ...
        fname,'fontsize',fsize,'fontweight',fontw,'backgroundcolor','w',varargin{:});
    % 

if loc(i)==3 || loc(i)==4
    set(htext(i),'horizontalalignment','r');
end

    set(hnewax(i),'visible','off','units',uni,'hittest','off');
    setappdata(hnewax(i),'NonDataObject',[]); % Used by DATACHILDREN.M.
    % This makes it
    % invulnerable to zooming etc
end;

return;
