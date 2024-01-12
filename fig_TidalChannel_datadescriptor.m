% -----------------------------
% A script to produce the Haro Starit Tidal Channel, ATOMIX benchmark 
% descriptor figure
% Ilker Fer, University of Bergen, Norway
% PI: Rolf Lueck, Rockland Scientific, Inc.
% -----------------------------

clear
close all
% add the folder with the dependencies to your path
addpath m_share
% point to where the benchmark data file is, and load
%data_path = 'C:\users\ngfif\Dropbox\ShearProbes\Data\VMP250_TidalChannel\' % ==> EDIT THIS
data_path = 'C:\ILKER\Dropbox\ShearProbes\Data\VMP250_TidalChannel\' % ==> EDIT THIS
%data_path = './'; % ==> EDIT THIS

save_plot_flag = 1;  % ==> Set to 1 if you want to save the figure as a PDF file
file_nc='VMP250_TidalChannel_024.nc';
plot_out_name = 'VMP250_TidalChannel.pdf'; % ==> name of the PDF file to save
% file_nc='VMP250_TidalChannel_024_cs.nc';
% plot_out_name = 'VMP250_TidalChannel_cs.pdf'; % ==> name of the PDF file to save

D = ATOMIX_load([data_path file_nc]); use D; clear D

%%
% convert time to seconds elapsed
TIME_ELAPSED_SEC_SLOW = (L1.TIME_SLOW - L1.TIME_SLOW(1))*(24*3600);
L1.TIME_ELAPSED_SEC   = (L1.TIME      - L1.TIME     (1))*(24*3600);
TIME_ELAPSED_SEC_EPSI = (L4.TIME      - L1.TIME     (1))*(24*3600);

SECTION_NUMBER_SLOW = interp1(L2.TIME, single(L2.SECTION_NUMBER),L1.TIME_SLOW);
PSPD_REL_SLOW       = interp1(L2.TIME, L2.PSPD_REL,L1.TIME_SLOW);
sec_slow            = find(SECTION_NUMBER_SLOW==1);
sec_fast            = find(L2.SECTION_NUMBER==1);

% some indices to bad (flag>0) data
ix_bad1     = find(L4.EPSI_FLAGS(:,1)>0);
ix_bad2     = find(L4.EPSI_FLAGS(:,2)>0);
ix_bad1_2     = find(L4.EPSI_FLAGS(:,1)+L4.EPSI_FLAGS(:,2)>0); % need these for scatter plot

% select estimate picks to show spectra
% good spectra, low eps
%pick1 = find(L4.FOM(:,1)<1.15 & L4.FOM(:,2)<1.15 & (L4.EPSI_FINAL<1e-7) );
% good spectra, med eps
%pick2 = find(L4.FOM(:,1)<1.15 & L4.FOM(:,2)<1.15 & (L4.EPSI_FINAL>1e-6 & L4.EPSI_FINAL<3.3e-6) );
% good spectra, hi eps
%pick3 = find(L4.FOM(:,1)<1.15 & L4.FOM(:,2)<1.15 & (L4.EPSI_FINAL>1e-5 & L4.EPSI_FINAL<3.3e-5) );
picks = [1 7]

%%
li = lines(6); % default colors of matlab... I'll force them just in case.
hyla=[]; hxla=[];

init_fig(14,17) % (WIDTH, HEIGHT) in cm, and outlay for the figure with subaxes
subaxis(7,1,1,'SV',0.01,'ML',.12,'MR',.1,'MB',.1,'MT',0.03);
% plot pressure time series on top panel
plot(TIME_ELAPSED_SEC_SLOW,L1.PRES_SLOW)
%hyla(1)=ylabel('$P\ \ [\mathrm{dbar}]$','interpreter', 'latex');
hyla(1)=ylabel({'$P$','$[\mathrm{dbar}]$'},'interpreter', 'latex');
% add W to axes 1
yyaxis right
plot(TIME_ELAPSED_SEC_SLOW,PSPD_REL_SLOW)
text([TIME_ELAPSED_SEC_SLOW(sec_slow(1)) TIME_ELAPSED_SEC_SLOW(sec_slow(end))],...
    [1.003 1.003],'\downarrow','fontsize',14,'fontw','b','horizontalal','center','verticalal','bottom')
ylabel('$W\ \ [\mathrm{m\,s^{-1}}]$','interpreter', 'latex')
set(gca,'ylim',[0 1]);

% plot shear probe time series with a given offset
subaxis(2);
offset = 10;
plot(L1.TIME_ELAPSED_SEC(sec_fast),L2.SHEAR(sec_fast,1),'color',li(1,:))
hold on
plot(L1.TIME_ELAPSED_SEC(sec_fast),L2.SHEAR(sec_fast,2)+offset,'color',li(2,:))
set(gca,'ylim',[-5 15]);
hyla(end+1)=ylabel({'$\partial u_i/\partial z$','$[\mathrm{s^{-1}}]$'},'interpreter', 'latex');

% plot acceleration or vibration time series
subaxis(3);
offset = 500;
plot(L1.TIME_ELAPSED_SEC(sec_fast),L2.VIB(sec_fast,1),'color',li(1,:))
hold on
plot(L1.TIME_ELAPSED_SEC(sec_fast),L2.VIB(sec_fast,2)+offset,'color',li(2,:))
set(gca,'ylim',[-600 800],'ytick',[-500:500:500]);
hyla(end+1)=ylabel({'$A_x,\ A_y$','$[\mathrm{counts}]$'},'interpreter', 'latex');

% Dissipation estimates from all probes and the final estimate
subaxis(4);
h1=semilogy(TIME_ELAPSED_SEC_EPSI,L4.EPSI(:,1),'o','color',li(1,:),'markerfacecolor', li(1,:), 'markersize',3);
hold on
h2=semilogy(TIME_ELAPSED_SEC_EPSI,L4.EPSI(:,2),'^','color',li(2,:),'markerfacecolor', li(2,:), 'markersize',3);
h3=semilogy(TIME_ELAPSED_SEC_EPSI,L4.EPSI_FINAL,'ks','markerfacecolor','w','markersize',3)
% % add red crosses for the flagged values
% semilogy(TIME_ELAPSED_SEC_EPSI(ix_bad1),L4.EPSI(ix_bad1,1),'rx','markersize',3);
% semilogy(TIME_ELAPSED_SEC_EPSI(ix_bad2),L4.EPSI(ix_bad2,2),'rx','markersize',3);


% mark the picked epsifinal for spectra 
%hvl = vlines([TIME_ELAPSED_SEC_EPSI(picks)],'k')
x=[TIME_ELAPSED_SEC_EPSI(picks)];   x=x(:)';  
ax=axis; y1=ax(3)*ones(size(x)); y2=ax(4)*ones(size(x));
x=[x;x]; y=[y1;y2];
hvl=plot(x,y); clear x y y1 y2 ax
uistack(hvl,'bottom');set(hvl,'color',[.5 .5 .5],'linew',1.2)

hyla(end+1)=ylabel({'$\varepsilon$','$[\mathrm{W\,kg^{-1}}]$'},'Interpreter','latex');
hleg=legend([h1(1) h2(1) h3(1)],{'$\varepsilon_1$','$\varepsilon_2$','$\varepsilon$'},...
    'location','southwest','fontsize',10,'Interpreter','latex');
%hleg.Box='off';
hleg.Orientation = 'horizontal';
hleg.ItemTokenSize(1)=12;
xlabel('$t\ \ [\mathrm{s}]$','interpreter', 'latex')

% tidy up time series axes
ax = flipud(findobj(gcf,'type','axes'));
set(ax(1:3),'XTickLabel',[]);
grid(ax,'on'); hold(ax,'on')
set(ax,'xlim',[0 200])

% set eps axis limits
ax(4).YLim =[.5e-7 2e-4];
ax(4).YTick=[1e-7 1e-6 1e-5 1e-4];
ax(4).YMinorGrid = 'off';ax(4).YMinorTick = 'off';


align_ylabel(hyla,-0.07)



%% add e1 vs e2 scatter plot
xli = [0.5e-7 1e-4]; % set limits after first test....
ax(5)=axes('position',[ax(1).Position(1) .06 .34 .37]);
ax(5).XLim =xli; 


h1=loglog(L4.EPSI(:,1),L4.EPSI(:,2),'^','color',li(2,:),'markerfacecolor',li(2,:),'markersize',3);
hold on
semilogy(L4.EPSI(ix_bad1_2,1),L4.EPSI(ix_bad1_2,2),'wx','markersize',2.5);
xli = get(gca,'xlim');
loglog(xli,xli,'w'); % 1-1 agreement
% bounds for stats uncertainty
diss_ratio_limit = 2.772; % this is ln(e1/e2) or ln(e1)-ln(2) times EPSI_STD
% EPSI_STD changes... just use the mean...
mean_e_std = mean(L4.EPSI_STD(:));
factor_1 = exp(diss_ratio_limit*mean_e_std);
patchx = [xli fliplr(xli) xli(1)];
patchy = [xli*(1/factor_1) fliplr(xli)*factor_1 xli(1)*(1/factor_1)];
hpa = patch(patchx,patchy,rgb('lightgray'),'LineStyle','none');
set(gca,'xlim',xli);
ch = get(gca,'children'); % this has the order
set(gca,'children',ch([2:length(ch) 1])); % put the patch back
axis square
ax(5).YLim =xli; ax(5).XLim=ax(5).YLim;
ax(5).XTick=[1e-7 1e-6 1e-5 1e-4];
ax(5).YTick=ax(5).XTick;
ax(5).XMinorGrid = 'off';ax(5).YMinorGrid = 'off';
grid(ax(5),'on')

xlabel('$\varepsilon_1\ \ [\mathrm{W\,kg^{-1}}]$','Interpreter','latex');
hyla(end+1)=ylabel('$\varepsilon_2\ \ [\mathrm{W\,kg^{-1}}]$','Interpreter','latex');

%% add shear spectra plot

ax(6)=axes('position',[0.56 ax(5).Position(2:4)]);

for I=1:2
    pick = picks(I);
    h1(1)=loglog(L3.KCYC(pick,:),L3.SH_SPEC_CLEAN(pick,:,1),'color',li(1,:));
    hold on
    h2(1)=loglog(L3.KCYC(pick,:),L3.SH_SPEC_CLEAN(pick,:,2),'color',li(2,:));
    ee = L4.EPSI_FINAL(pick);
    [nas,~] = lueck_spectrum(ee, L4.KVISC(pick),L3.KCYC(pick,:));
    h3(1)=loglog(L3.KCYC(pick,:),nas,'color',[.6 .6 .6],'linew',1.2);
end

set(gca,'xlim',[1 150],'ylim',[1e-6 2e-1])
grid
xlabel('$k\ \ [\mathrm{cpm}]$','Interpreter','latex');
ylabel('$\Psi\ [\mathrm{s^{-2}\,cpm^{-1}}]$','Interpreter','latex');

hleg1 = legend([h1 h2 h3],{'$\Psi_{1c}$','$\Psi_{2c}$','Lueck'},...
    'location','southwest','Interpreter','latex');
hleg1.ItemTokenSize(1)=20;
axis square
ax(6).XMinorGrid = 'off';ax(6).YMinorGrid = 'off';




%% some cosmetics; putting a, b, etc.
finish_fig(ax,[],[],[],[],[],'times');

if save_plot_flag
    exportgraphics(gcf,plot_out_name,'ContentType', 'vector');
end
%%
