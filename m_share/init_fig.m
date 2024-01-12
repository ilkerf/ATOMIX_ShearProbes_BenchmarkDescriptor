function init_fig(WIDTH,HEIGHT);
% init_fig(WIDTH,HEIGHT);
% SEE ALSO finish_fig.m

if nargin<1
    WIDTH=12; HEIGHT=8;
end

d = 0.2;
figure('units','centimeters','color','w',...
    'position',[d d WIDTH-d HEIGHT-d],...
    'Paperposition',[d d WIDTH-d HEIGHT-d],...
    'PaperSize',[WIDTH+2*d HEIGHT+2*d]);


