global pa0 cm0 cbio0131 cbio01 cbio0141 cbio0132 cbio02 cbio0142 ...
    cbio0133 cbio03 cbio0143 cm013 pa013 cm014 pa014 C14prog ...
    akba1 akba2 akba3 ...
    akab1 akab2 akab3 epsi akam ak 

% Script showing examples of plotting model simulation output
% Heather Graven 2020

% Load model output
load('SSPsims_200724084931.mat')
C14prog=0; % change to zero so that BDprecalc does not attempt to load ex14sourcef input file

for t=1:length(param(:,1))
    % specify parameters that were used for each model run and calculate
    % reference values for each model run
    ak=param(t,1); akam=1./param(t,2); epsi=param(t,3);
    akab1=1/param(t,4); akab2=1/param(t,5); akab3=1/param(t,6);
    akba1=1/param(t,7); akba2=1/param(t,8); akba3=1/param(t,9);
    BDprecalc;
    cbio01_p(1,t)=cbio01; cbio0131_p(1,t)=cbio0131; cbio0141_p(1,t)=cbio0141;
    cbio02_p(1,t)=cbio02; cbio0132_p(1,t)=cbio0132; cbio0142_p(1,t)=cbio0142;
    cbio03_p(1,t)=cbio03; cbio0133_p(1,t)=cbio0133; cbio0143_p(1,t)=cbio0143;
    cm0_p(1,t)=cm0; cm013_p(1,t)=cm013; cm014_p(1,t)=cm014;
    pa0_p(1,t)=pa0; pa013_p(1,t)=pa013; pa014_p(1,t)=pa014;
end

%%%%%% calculate simulated d13C and D14C in each pool
% Atmosphere
d14C(:,1,:,:)=((squeeze(SSPY(:,1+(3-1)*49,:,:))+repmat(pa014_p,[601 1 6]))./(squeeze(SSPY(:,1,:,:))+repmat(pa0_p,[601 1 6]))-1)*1000;
d13C(:,1,:,:)=((squeeze(SSPY(:,1+(2-1)*49,:,:))+repmat(pa013_p,[601 1 6]))./(squeeze(SSPY(:,1,:,:))+repmat(pa0_p,[601 1 6]))/rs-1)*1000;
D14C(:,1,:,:)=d14C(:,1,:,:)-2*(d13C(:,1,:,:)+25).*(1+d14C(:,1,:,:)/1000);
% Ocean
for p=2:44
    d14C(:,p,:,:)=((squeeze(SSPY(:,p+(3-1)*49,:,:))+repmat(cm014_p,[601 1 6]))./(squeeze(SSPY(:,p,:,:))+repmat(cm0_p,[601 1 6]))-1)*1000;
    d13C(:,p,:,:)=((squeeze(SSPY(:,p+(2-1)*49,:,:))+repmat(cm013_p,[601 1 6]))./(squeeze(SSPY(:,p,:,:))+repmat(cm0_p,[601 1 6]))/rs-1)*1000;
    D14C(:,p,:,:)=d14C(:,p,:,:)-2*(d13C(:,p,:,:)+25).*(1+d14C(:,p,:,:)/1000);
end
% Biosphere pool 1
d14C(:,45,:,:)=((squeeze(SSPY(:,45+(3-1)*49,:,:))+repmat(cbio0141_p,[601 1 6]))./(squeeze(SSPY(:,45,:,:))+repmat(cbio01_p,[601 1 6]))-1)*1000;
d13C(:,45,:,:)=((squeeze(SSPY(:,45+(2-1)*49,:,:))+repmat(cbio0131_p,[601 1 6]))./(squeeze(SSPY(:,45,:,:))+repmat(cbio01_p,[601 1 6]))/rs-1)*1000;
D14C(:,45,:,:)=d14C(:,45,:,:)-2*(d13C(:,45,:,:)+25).*(1+d14C(:,45,:,:)/1000);
% Biosphere pool 2
d14C(:,46,:,:)=((squeeze(SSPY(:,46+(3-1)*49,:,:))+repmat(cbio0142_p,[601 1 6]))./(squeeze(SSPY(:,46,:,:))+repmat(cbio02_p,[601 1 6]))-1)*1000;
d13C(:,46,:,:)=((squeeze(SSPY(:,46+(2-1)*49,:,:))+repmat(cbio0132_p,[601 1 6]))./(squeeze(SSPY(:,46,:,:))+repmat(cbio02_p,[601 1 6]))/rs-1)*1000;
D14C(:,46,:,:)=d14C(:,46,:,:)-2*(d13C(:,46,:,:)+25).*(1+d14C(:,46,:,:)/1000);
% Biosphere pool 3
d14C(:,47,:,:)=((squeeze(SSPY(:,47+(3-1)*49,:,:))+repmat(cbio0143_p,[601 1 6]))./(squeeze(SSPY(:,47,:,:))+repmat(cbio03_p,[601 1 6]))-1)*1000;
d13C(:,47,:,:)=((squeeze(SSPY(:,47+(2-1)*49,:,:))+repmat(cbio0133_p,[601 1 6]))./(squeeze(SSPY(:,47,:,:))+repmat(cbio03_p,[601 1 6]))/rs-1)*1000;
D14C(:,47,:,:)=d14C(:,47,:,:)-2*(d13C(:,47,:,:)+25).*(1+d14C(:,47,:,:)/1000);

% max, min, mid values
maxD14C=squeeze(max(D14C,[],3));
minD14C=squeeze(min(D14C,[],3));
midD14C=(maxD14C+minD14C)/2;

maxd13C=squeeze(max(d13C,[],3));
mind13C=squeeze(min(d13C,[],3));
midd13C=(maxd13C+mind13C)/2;

% Define colors for plotting 6 SSPs
sspco=[[25 76 156]/255;
    [62 153 17]/255;
    [.45 .26 0.26];
    [153 153 0]/255;
    [.6 .2 0];
    [.85 .33 0.1]];

% plot atmospheric D14C
figure; hold on
for rc=1:6
    fill([2005:2100 fliplr(2005:2100)],[squeeze(minD14C(506:end,1,rc))' ...
        fliplr(squeeze(maxD14C(506:end,1,rc))')],sspco(rc,:),'LineStyle','none')
end
plot(inputdata.C14time,inputdata.C14data,'k-','LineWidth',1.5);
set(gca,'Fontsize',11)
xlim([1940 2100])
ylim([-300 800])
legend('SSP1-1.9','SSP1-2.6','SSP2-4.5','SSP3-7.0','SSP5-3.4os','SSP5-8.5','Observed')
plot([1940 2100],[0 0],'k-')
ylabel('\Delta^{14}CO_2')
ax1=gca;
ax2 = axes('Position',get(ax1,'Position'),...
    'XAxisLocation','top','XTickLabel',[],...
    'YAxisLocation','right',...
    'Color','none',...
    'XColor','k','YColor','k');
set(ax2,'YTick',(exp((-3000:500:0)/8033)-1)*1000,'YTickLabel',{'3000';'2500';'2000';'1500';'1000';'500';'0'})
xlim([1940 2100])
ylim([-300 800])

% plot atmospheric d13C
figure; hold on
for rc=1:6
    fill([2005:2100 fliplr(2005:2100)],[squeeze(mind13C(506:end,1,rc))' ...
        fliplr(squeeze(maxd13C(506:end,1,rc))')],sspco(rc,:),'LineStyle','none')
end
fill([1500:2005 fliplr(1500:2005)],[squeeze(mind13C(1:506,1,rc))' ...
        fliplr(squeeze(maxd13C(1:506,1,rc))')],[.7 .7 .7],'LineStyle','none')
plot(inputdata.C13time,inputdata.C13data,'k-','LineWidth',1.5);
set(gca,'Fontsize',11,'layer','top','box','on')
xlim([1940 2100])
ylim([-14 -6])
ylabel('\delta^{13}CO_2')

% create tables of atmospheric D14C and d13C outputs
etb=(2005:2100)'; detb=(2005:2100)';
for rc=1:6
    etb=[etb squeeze(midD14C(506:end,1,rc)) squeeze(minD14C(506:end,1,rc)) squeeze(maxD14C(506:end,1,rc))];
    detb=[detb squeeze(midd13C(506:end,1,rc)) squeeze(mind13C(506:end,1,rc)) squeeze(maxd13C(506:end,1,rc))];
end

% Additional plots
% Plot D14C and d13C in other pools
for r=1:6
    figure; hold on
    plot(1500:2100,squeeze(midD14C(:,1,r)),'k-','LineWidth',1.5);
    
    fill([1500:2100 fliplr(1500:2100)],[squeeze(minD14C(:,45,r))' fliplr(squeeze(maxD14C(:,45,r))')],'g','LineStyle','none')
    fill([1500:2100 fliplr(1500:2100)],[squeeze(minD14C(:,2,r))' fliplr(squeeze(maxD14C(:,2,r))')],'c','LineStyle','none')
    fill([1500:2100 fliplr(1500:2100)],[squeeze(minD14C(:,13,r))' fliplr(squeeze(maxD14C(:,13,r))')],'b','LineStyle','none')
    fill([1500:2100 fliplr(1500:2100)],[squeeze(minD14C(:,25,r))' fliplr(squeeze(maxD14C(:,25,r))')],'m','LineStyle','none')
    fill([1500:2100 fliplr(1500:2100)],[squeeze(minD14C(:,44,r))' fliplr(squeeze(maxD14C(:,44,r))')],[.5 .5 .5],'LineStyle','none')
    
    xlim([1940 2100])
    legend('atmosphere','bio box 1','mixed layer','ocean 300m','ocean 600m','ocean 3500m')
    title(['SSP' SSP{r}])
    ylim([-300 800])
    plot([1940 2100],[0 0],'k-')
    set(gca,'XTick',1940:20:2100,'Fontsize',11,'box','on','layer','top')
    ylabel('\Delta^{14}C')
end

for r=1:6
    figure; hold on
    plot(1500:2100,squeeze(midd13C(:,1,r)),'k-','LineWidth',1.5);
    
    fill([1500:2100 fliplr(1500:2100)],[squeeze(mind13C(:,45,r))' fliplr(squeeze(maxd13C(:,45,r))')],'g','LineStyle','none')
    fill([1500:2100 fliplr(1500:2100)],[squeeze(mind13C(:,2,r))' fliplr(squeeze(maxd13C(:,2,r))')],'c','LineStyle','none')
    fill([1500:2100 fliplr(1500:2100)],[squeeze(mind13C(:,13,r))' fliplr(squeeze(maxd13C(:,13,r))')],'b','LineStyle','none')
    fill([1500:2100 fliplr(1500:2100)],[squeeze(mind13C(:,25,r))' fliplr(squeeze(maxd13C(:,25,r))')],'m','LineStyle','none')
    fill([1500:2100 fliplr(1500:2100)],[squeeze(mind13C(:,44,r))' fliplr(squeeze(maxd13C(:,44,r))')],[.5 .5 .5],'LineStyle','none')
    
    xlim([1940 2100])
    legend('atmosphere','bio box 1','mixed layer','ocean 300m','ocean 600m','ocean 3500m')
    title(['SSP' SSP{r}])
    plot([1940 2100],[0 0],'k-')
    set(gca,'XTick',1940:20:2100,'Fontsize',11,'box','on','layer','top')
    ylabel('\delta^{13}C')
end

%calculate total atoms C14, in units 10^26 atoms
atoms14=zeros([601 45 length(param(:,1)) 6]);
atoms14(:,1,:,:)=(squeeze(SSPY(:,1+(3-1)*49,:,:))+repmat(pa014_p,[601 1 6]))/gtc2ppm*1e15/12*R14s*6.02e-3;
atoms14(:,2,:,:)=(squeeze(SSPY(:,2+(3-1)*49,:,:))+repmat(cm014_p,[601 1 6]))*R14s*hm*aoc/1000*6.02e-3;
for i=3:39; atoms14(:,i,:,:)=(squeeze(SSPY(:,i+(3-1)*49,:,:))+repmat(cm014_p,[601 1 6]))*R14s*h1*aoc/1000*6.02e-3; end
for i=40:44; atoms14(:,i,:,:)=(squeeze(SSPY(:,i+(3-1)*49,:,:))+repmat(cm014_p,[601 1 6]))*R14s*h2*aoc/1000*6.02e-3; end
atoms14(:,45,:,:)=(squeeze(SSPY(:,45+(3-1)*49,:,:))+repmat(cbio0141_p,[601 1 6]))*R14s*1e15/12*6.02e-3;
atoms14(:,46,:,:)=(squeeze(SSPY(:,46+(3-1)*49,:,:))+repmat(cbio0142_p,[601 1 6]))*R14s*1e15/12*6.02e-3;
atoms14(:,47,:,:)=(squeeze(SSPY(:,47+(3-1)*49,:,:))+repmat(cbio0143_p,[601 1 6]))*R14s*1e15/12*6.02e-3;
atoms14oc=squeeze(sum(atoms14(:,2:44,:,:),2));

% plot 14C inventories in atmosphere, ocean and biosphere
figure('Units','Inches','Position',[1    1   12.8021    3.2813]);
subplot(1,3,1);
title('Atmosphere')
hold on
for rc=1:6
    fill([1500:2100 fliplr(1500:2100)],[squeeze(max(atoms14(:,1,:,rc)-...
        repmat(atoms14(451,1,:,rc),[601 1 1 1]),[],3))' ...
        fliplr(squeeze(min(atoms14(:,1,:,rc)-repmat(atoms14(451,1,:,rc),[601 1 1 1])...
        ,[],3))')],sspco(rc,:),'LineStyle','none')
end
plot(1500:2005,squeeze(atoms14(1:506,1,1,1))-squeeze(atoms14(451,1,1,1)),'k-','LineWidth',1)
plot([1940 2100],[0 0],'k-')
xlim([1945 2100])
ylim([-100 600])
ylabel('14C inventory')
set(gca,'box','on','XTick',1950:25:2100)
% Ocean
subplot(1,3,2);
title('Ocean')
hold on
for rc=1:6
    fill([1500:2100 fliplr(1500:2100)],[squeeze(max(atoms14oc(:,:,rc)-...
        repmat(atoms14oc(451,:,rc),[601 1 1]),[],2))' ...
        fliplr(squeeze(min(atoms14oc(:,:,rc)-repmat(atoms14oc(451,:,rc),[601 1 1]),...
        [],2))')],sspco(rc,:),'LineStyle','none')
end
plot(1500:2005,squeeze(atoms14(1:506,1,1,1))-squeeze(atoms14(451,1,1,1)),'k-','LineWidth',1)
plot([1940 2100],[0 0],'k-')
xlim([1945 2100])
ylim([-100 600])
ylabel('14C inventory')
set(gca,'box','on','XTick',1950:25:2100)
% Biosphere
subplot(1,3,3);
title('Biosphere')
hold on
for rc=1:6
    fill([1500:2100 fliplr(1500:2100)],[squeeze(max(sum(atoms14(:,45:47,:,rc),2)-...
        repmat(sum(atoms14(451,45:47,:,rc),2),[601 1 1 1]),[],3))' ...
        fliplr(squeeze(min(sum(atoms14(:,45:47,:,rc),2)-repmat(sum(atoms14(451,45:47,:,rc),2),[601 1 1 1]),...
        [],3))')],sspco(rc,:),'LineStyle','none')
end
plot(1500:2005,squeeze(atoms14(1:506,1,1,1))-squeeze(atoms14(451,1,1,1)),'k-','LineWidth',1)
plot([1940 2100],[0 0],'k-')
xlim([1945 2100])
ylim([-100 600])
ylabel('14C inventory')
set(gca,'box','on','XTick',1950:25:2100)
set(gcf,'PaperPositionMode','auto','PaperSize',[12.8021    3.2813]);

% Plot total 14C in carbon cycle
figure; hold on
plot(1500:2100,squeeze(sum(atoms14(:,1:47,:,rc),2)-repmat(sum(atoms14(451,1:47,:,rc),2),[601 1 1 1])))
plot([1990 2000 2000 1990 1990],615+35*[-1 -1 1 1 -1],'k-')
xlim([1940 2100])

% From Naegler et al. 2009, ocean bomb C14 inventory of 374 ± 98 ? 10^26 atoms
% 1995 minus 1950
figure; hold on
plot(1500:2100,squeeze(atoms14oc(:,:,1)-repmat(atoms14oc(451,:,1),[601 1 1])))
plot([1994 1996 1996 1994 1994],374+98*[-1 -1 1 1 -1],'k-')
xlim([1940 2005])

%calculate total CO2 in ocean
TotOcAnthC=zeros([601 length(param(:,1)) 6]); 
for r=1:6
    TotOcAnthC(:,:,r)=aoc/1000*12/1e15*(hm*(squeeze(SSPY(:,2,:,r))+repmat(cm0_p,[601 1])) ...
        +squeeze(sum(h1*(squeeze(SSPY(:,3:39,:,r))+squeeze(repmat(reshape(cm0_p,[1 1 length(cm0_p)]),[601 37 1]))),2)) ...
        +squeeze(sum(h2*(squeeze(SSPY(:,40:44,:,r))+squeeze(repmat(reshape(cm0_p,[1 1 length(cm0_p)]),[601 5 1]))),2)));
end
% From Sabine et al. 2004, ocean co2 inventory for 1800 to 1994 of 118 +/- 19 petagrams of carbon
% appears too low because simulation started in 1850
figure; hold on
plot(1500:2100,squeeze(TotOcAnthC(:,:,1)-repmat(TotOcAnthC(301,:,1),[601 1 1])))
plot([1993 1995 1995 1993 1993],118+19*[-1 -1 1 1 -1],'k-')
xlim([1940 2005])

% plot land carbon stocks
for rc=1:6
    figure
    subplot(1,4,1); hold on
    plot(1500:2100,squeeze(sum(SSPY(:,45,:,rc),2))+repmat(cbio01_p(:)',[601 1]))
    xlim([1850 2100])
    title('Box 1')
    subplot(1,4,2); hold on
    plot(1500:2100,squeeze(sum(SSPY(:,46,:,rc),2))+repmat(cbio02_p(:)',[601 1]))
    xlim([1850 2100])
    title('Box 2')
    subplot(1,4,3); hold on
    plot(1500:2100,squeeze(sum(SSPY(:,47,:,rc),2))+repmat(cbio03_p(:)',[601 1]))
    xlim([1850 2100])
    title('Box 3')
    subplot(1,4,4); hold on
    plot(1500:2100,squeeze(sum(SSPY(:,45:47,:,rc),2)))
    xlim([1850 2100])
    title(['Total C anomaly, SSP' SSP{rc}])
end


