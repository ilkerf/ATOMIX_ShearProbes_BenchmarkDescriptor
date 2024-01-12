function hla=finish_fig(ax,varargin);
% hla=finish_fig(ax,hlabs,sublabloc,fonts,lab_fonts,ticklen,fontn,fontw);
% SEE ALSO init_fig.m
% ax is required, others have defaults...
% fonts=9; lab_fonts=10; ticklen=[.02 0.01]; font is arial
% If you don't want subaxis labels use 0 for sublablac
% fontw = 'b' is for subaxis label
% sublabloc are  1  top left (default), 2  bottom left, 3  bottom right, 4  top right
% for 4 axes, can assign different locations: e.g. [1 2 1 1]


% set defaults for optional inputs
optargs = {[] ones(size(ax)) 9 10 [.02 0.01] 'arial' 'b'};
% replace defaults by inputs
% account for eventual empty inputs to activate default
newVals = cellfun(@(x) ~isempty(x), varargin);
optargs(newVals) = varargin(newVals);
[hlabs,sublabloc,fonts,lab_fonts,ticklen,fontn,fontw] = optargs{:};
%---------------------------------



set(ax,'fontsize',fonts,'fontname',fontn,'ticklen',ticklen);
set(hlabs,'fontname',fontn,'fontsize',lab_fonts);

if length(ax)>1 && sublabloc(1)>0
hla=sublabel(ax,[],[],sublabloc,'fontname',fontn,'fontsize',10,'fontweight',fontw,'margin',.1);
else
    hla=[];
end
