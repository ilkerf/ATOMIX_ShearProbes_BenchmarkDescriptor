% -----------------------------
% A script to produce the Nemo_MR1000_Minas_Passage_InStream, ATOMIX benchmark 
% descriptor figure
% code template: Ilker Fer, University of Bergen, Norway
% PI: Rolf Lueck, Rockland Scientific, Inc.
% -----------------------------

clear
%close all
% add the folder with the dependencies to your path
%addpath m_share
% point to where the benchmark data file is, and load
data_path = 'C:\ILKER\Dropbox\ShearProbes\Data\Nemo_MR1000_Minas_Passage\' ;% ==> EDIT THIS
%data_path = './'; % ==> EDIT THIS

save_plot_flag = 1; % 0;  % ==> Set to 1 if you want to save the figure as a PDF file
file_nc='Nemo_MR1000_Minas_Passage_InStream.nc';
plot_out_name = 'Nemo_MR1000_Minas_Passage.pdf'; % ==> name of the PDF file to save

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
ix_bad3     = find(L4.EPSI_FLAGS(:,3)>0);
ix_bad4     = find(L4.EPSI_FLAGS(:,4)>0);
ix_bad1_2     = find(L4.EPSI_FLAGS(:,1)+L4.EPSI_FLAGS(:,2)>0); % need these for scatter plot
ix_bad1_3     = find(L4.EPSI_FLAGS(:,1)+L4.EPSI_FLAGS(:,3)>0);
ix_bad1_4     = find(L4.EPSI_FLAGS(:,1)+L4.EPSI_FLAGS(:,4)>0);
ix_allbad   = find(sum(L4.EPSI_FLAGS,2)>0);



% select estimate picks to show spectra
picks = 28; % show only 1 pick, we have 4 probes!

%%
li = lines(6); % default colors of matlab... I'll force them just in case.
hyla=[]; hxla=[];

init_fig(14,17) % (WIDTH, HEIGHT) in cm, and outlay for the figure with subaxes
subaxis(7,1,1,'SV',0.01,'ML',.12,'MR',.1,'MB',.1,'MT',0.03);
% plot pressure time series on top panel
plot(TIME_ELAPSED_SEC_SLOW,L1.PRES_SLOW); 
set(gca,'ylim',[42 44],'ytick',[42 43 44]);
hyla(1)=ylabel({'$P$','$[\mathrm{dbar}]$'},'interpreter', 'latex');
% add W to axes 1
yyaxis right
plot(TIME_ELAPSED_SEC_SLOW,PSPD_REL_SLOW)
text([TIME_ELAPSED_SEC_SLOW(sec_slow(1)) TIME_ELAPSED_SEC_SLOW(sec_slow(end))],...
    [2.003 2.003],'\downarrow','fontsize',14,'fontw','b','horizontalal','center','verticalal','bottom')
ylabel('$U\ \ [\mathrm{m\,s^{-1}}]$','interpreter', 'latex')
set(gca,'ylim',[0 2]);

% plot shear probe time series with a given offset
subaxis(2);
offset = 50;
plot(L1.TIME_ELAPSED_SEC(sec_fast),L2.SHEAR(sec_fast,1),'color',li(1,:))
hold on
plot(L1.TIME_ELAPSED_SEC(sec_fast),L2.SHEAR(sec_fast,2)+offset,'color',li(2,:))
plot(L1.TIME_ELAPSED_SEC(sec_fast),L2.SHEAR(sec_fast,3)+2*offset,'color',li(3,:))
plot(L1.TIME_ELAPSED_SEC(sec_fast),L2.SHEAR(sec_fast,4)+3*offset,'color',li(4,:))
set(gca,'ylim',[-50 210],'ytick',[0 100 200]);
hyla(end+1)=ylabel({'$\partial u_i/\partial x$','$[\mathrm{s^{-1}}]$'},'interpreter', 'latex');

% plot acceleration or vibration time series
subaxis(3);
offset = 4000;
plot(L1.TIME_ELAPSED_SEC(sec_fast),L2.VIB(sec_fast,1),'color',li(1,:))
hold on
plot(L1.TIME_ELAPSED_SEC(sec_fast),L2.VIB(sec_fast,2)+offset,'color',li(2,:))
set(gca,'ylim',[-5000 6000],'ytick',[-2000:2000:4000]);
hyla(end+1)=ylabel({'$A_z,\ A_y$','$[\mathrm{counts}]$'},'interpreter', 'latex');

% Dissipation estimates from all probes and the final estimate
subaxis(4);
h1=semilogy(TIME_ELAPSED_SEC_EPSI,L4.EPSI(:,1),'o','color',li(1,:),'markerfacecolor', li(1,:), 'markersize',3);
hold on
h2=semilogy(TIME_ELAPSED_SEC_EPSI,L4.EPSI(:,2),'^','color',li(2,:),'markerfacecolor', li(2,:), 'markersize',3);
h3=semilogy(TIME_ELAPSED_SEC_EPSI,L4.EPSI(:,3),'o','color',li(3,:),'markerfacecolor', li(3,:), 'markersize',3);
h4=semilogy(TIME_ELAPSED_SEC_EPSI,L4.EPSI(:,4),'^','color',li(4,:),'markerfacecolor', li(4,:), 'markersize',3);
h5=semilogy(TIME_ELAPSED_SEC_EPSI,L4.EPSI_FINAL,'ks','markerfacecolor','w','markersize',3);
% add red crosses for the flagged values
%semilogy(TIME_ELAPSED_SEC_EPSI(ix_bad1),L4.EPSI(ix_bad1,1),'rx','markersize',3);
%semilogy(TIME_ELAPSED_SEC_EPSI(ix_bad2),L4.EPSI(ix_bad2,2),'rx','markersize',3);
%semilogy(TIME_ELAPSED_SEC_EPSI(ix_bad3),L4.EPSI(ix_bad3,3),'rx','markersize',3);
%semilogy(TIME_ELAPSED_SEC_EPSI(ix_bad4),L4.EPSI(ix_bad4,4),'rx','markersize',3);

% mark the picked epsifinal for spectra 
set(gca,'ylim',[1e-5 3e-3]);
x=[TIME_ELAPSED_SEC_EPSI(picks)];   x=x(:)';  
ax=axis; y1=ax(3)*ones(size(x)); y2=ax(4)*ones(size(x));
x=[x;x]; y=[y1;y2];
hvl=plot(x,y); clear x y y1 y2 ax
uistack(hvl,'bottom');set(hvl,'color',[.5 .5 .5],'linew',1.2)



hyla(end+1)=ylabel({'$\varepsilon$','$[\mathrm{W\,kg^{-1}}]$'},'Interpreter','latex');
hleg=legend([h1(1) h2(1) h3(1) h4(1) h5(1)],...
    {'$\varepsilon_1$','$\varepsilon_2$','$\varepsilon_3$','$\varepsilon_4$','$\varepsilon$'},...
    'location','southeast','fontsize',10,'Interpreter','latex');
%hleg.Box='off';
hleg.Orientation = 'horizontal';
hleg.ItemTokenSize(1)=10;
xlabel('$t\ \ [\mathrm{s}]$','interpreter', 'latex')

% tidy up time series axes
ax = flipud(findobj(gcf,'type','axes'));
set(ax(1:3),'XTickLabel',[]);
grid(ax,'on'); hold(ax,'on')
set(ax,'xlim',[0 310])

% set eps axis limits
ax(4).YLim =[1e-5 3e-3];
ax(4).YTick=[1e-5 1e-4 1e-3];
ax(4).YMinorGrid = 'off';ax(4).YMinorTick = 'off';


align_ylabel(hyla,-0.07)



%% add e1 vs e2 scatter plot
xli = [1e-4 3e-3]; % set limits after first test....
ax(5)=axes('position',[ax(1).Position(1) .06 .34 .37]);
ax(5).XLim =xli; 


h1=loglog(L4.EPSI(:,1),L4.EPSI(:,2),'^','color',li(2,:),'markerfacecolor',li(2,:),'markersize',3);
hold on
h2=loglog(L4.EPSI(:,1),L4.EPSI(:,3),'o','color',li(3,:),'markerfacecolor',li(3,:),'markersize',3);
h3=loglog(L4.EPSI(:,1),L4.EPSI(:,4),'^','color',li(4,:),'markerfacecolor',li(4,:),'markersize',3);
semilogy(L4.EPSI(ix_bad1_2,1),L4.EPSI(ix_bad1_2,2),'wx','markersize',2.5);
semilogy(L4.EPSI(ix_bad1_3,1),L4.EPSI(ix_bad1_3,3),'wx','markersize',2.5);
semilogy(L4.EPSI(ix_bad1_4,1),L4.EPSI(ix_bad1_4,4),'wx','markersize',2.5);

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
%ax(5).XTick=[1e-5 1e-4 1e-3 1e-2];
ax(5).XTick=[1e-4 1e-3];
ax(5).YTick=ax(5).XTick;
ax(5).XMinorGrid = 'off';ax(5).YMinorGrid = 'off';
grid(ax(5),'on')

xlabel('$\varepsilon_1\ \ [\mathrm{W\,kg^{-1}}]$','Interpreter','latex');
hyla(end+1)=ylabel('$\varepsilon_2, \varepsilon_3, \varepsilon_4 \ [\mathrm{W\,kg^{-1}}]$','Interpreter','latex');
hleg=legend([h1(1) h2(1) h3(1)],...
    {'$(\varepsilon_1,\varepsilon_2)$','$(\varepsilon_1,\varepsilon_3)$','$(\varepsilon_1,\varepsilon_4)$'},...
    'location','southeast','fontsize',10,'Interpreter','latex');
hleg.ItemTokenSize(1)=20;
%% add shear spectra plot
clear h1 h2 h3 h4
ax(6)=axes('position',[0.56 ax(5).Position(2:4)]);

for I=1:length(picks)
    pick = picks(I);
    h1(1)=loglog(L3.KCYC(pick,:),L3.SH_SPEC_CLEAN(pick,:,1),'color',li(1,:));
    hold on
    h2(1)=loglog(L3.KCYC(pick,:),L3.SH_SPEC_CLEAN(pick,:,2),'color',li(2,:));
    h3(1)=loglog(L3.KCYC(pick,:),L3.SH_SPEC_CLEAN(pick,:,3),'color',li(3,:));
    h4(1)=loglog(L3.KCYC(pick,:),L3.SH_SPEC_CLEAN(pick,:,4),'color',li(4,:));
    ee = L4.EPSI_FINAL(pick);
    [nas,~] = lueck_spectrum(ee, L4.KVISC(pick),L3.KCYC(pick,:));
    h5(1)=loglog(L3.KCYC(pick,:),nas,'color',[.6 .6 .6],'linew',1.2);
end

set(gca,'xlim',[1 150],'ylim',[1e-2 1])
grid on
xlabel('$k\ \ [\mathrm{cpm}]$','Interpreter','latex');
ylabel('$\Psi\ [\mathrm{s^{-2}\,cpm^{-1}}]$','Interpreter','latex');

hleg1 = legend([h1 h2 h3 h4 h5],{'$\Psi_{1}$','$\Psi_{2}$','$\Psi_{3}$','$\Psi_{4}$','Lueck'},...
    'location','southeast','Interpreter','latex', 'NumColumns', 2);%, 'box', 'off');
hleg1.ItemTokenSize(1)=15;
axis square
ax(6).XMinorGrid = 'off';ax(6).YMinorGrid = 'off';




%% some cosmetics; putting a, b, etc.
%finish_fig(ax,[],[],[],[],[],'times');
finish_fig(ax,[],[2 1 2 2 1 1],[],[],[],'times')

if save_plot_flag
    exportgraphics(gcf,plot_out_name,'ContentType', 'vector');
end
%%
