% Generate analysis/plots for Temperature impacts paper
clearvars -except S;
close all; clc;

RELOAD=1;
if RELOAD==1
    load(append('data','/','S9.mat'))%load('S9.mat'); # Request the data files from the authors
    
end

cd RCP
[ti, Tmax, Tmean, Tmin]=temperatureProjectionsCMIP5();
cd ..

QC={S.QC};
ID=find(strcmp(QC,'Everything seems ok.'));

% remove some stations manually
ID(strcmp({S(ID).name}','THORNBROUGH'))=[];
ID(strcmp({S(ID).name}','Buckskin Joe'))=[];
ID(strcmp({S(ID).name}','Saint Elmo'))=[];
ID(strcmp({S(ID).name}','WAPELLO'))=[];
ID(strcmp({S(ID).name}','HYDEN 4 E'))=[];
ID(strcmp({S(ID).name}','ANDERSON 3W'))=[];
ID(strcmp({S(ID).name}','AVA 6NW'))=[];
ID(strcmp({S(ID).name}','JOPLIN 24 N'))=[];
ID(strcmp({S(ID).name}','PARKS SPRING RCH'))=[];
ID(strcmp({S(ID).name}','CHARLOTTE'))=[];
ID(strcmp({S(ID).name}','MT. PLEASANT'))=[];
ID(strcmp({S(ID).name}','Swift Creek'))=[];
ID(strcmp({S(ID).name}','LAC VIEUX DESERT'))=[];

ID([S(ID).yearsWithData]<=50)=[];
len=length(ID);

% Stations of interests 
cities={'USW00022521','Honolulu'
    'USW00023272','San Francisco'
    'USW00024233','Seattle'
    'USW00025309','Juneau'
    'USW00013960','Dallas'
    'USW00094846','Chicago'
    'USW00014739','Boston'
    'USW00013743','Washington D.C.'
    };
ID_cities=[];
for i =1:size(cities,1)
    ID_cities(end+1)=find(strcmp({S(ID).id},cities{i,1}));
end
ri=1/3;
% anonymous functions to work with GEV data
F=@(x,mu,sig,k) gevcdf(x,k,sig,mu); % Cumulative probablity distribution function
E=@(x,mu,sig,k) 1-F(x,mu,sig,k);    % exceedance probablity function
E_inc=@(x,mu,sig,k,SLR) E(x,mu+SLR,sig,k)./E(x,mu,sig,k); % factor of increase in exceedance of probability
x_inv=@(E,mu,sig,k) gevinv(1-E,k,sig,mu);                 % calculate event level x based on exceedace probability
Tr=@(x,mu,sig,k) ri./E(x,mu,sig,k);                       % return period
x_Tr=@(yr,mu,sig,k) gevinv(1-(ri./yr),k,sig,mu);          % calculate event level x based on return period
Tr_pct_red=@(TR,mu,sig,k,SLR) (TR-Tr(x_Tr(TR,mu,sig,k),mu+SLR,sig,k))./TR; % calculate the percent reduction in return period

DT=@(mu,sig,k) x_Tr(50,mu,sig,k)-x_Tr(1,mu,sig,k); % difference in temperature level between 50-year and 1-year event

DeltaT_50yr_1yr=NaN(len,1);
DeltaT_1yr_SUMTEMP=NaN(len,1);

for i=1:len
    DeltaT_50yr_1yr(i)=DT(S(ID(i)).mu,S(ID(i)).sigma,S(ID(i)).k);
    DeltaT_1yr_SUMTEMP(i)=x_Tr(1,S(ID(i)).mu,S(ID(i)).sigma,S(ID(i)).k)-S(ID(i)).summer_avg;
end

% % solve for sigma as a function of K and DWL
mu=40;
k=(-1:0.02:1)';
DTi=(0.01:0.1:50)';

[K,DTi]=meshgrid(k,DTi);

[m,n]=size(K);
SIG=NaN(m,n);

for i=1:m
    for j=1:n
        SIG(i,j)=fzero(@(sig) DT(mu,sig*(sig>0)+1e-6,K(i,j))-DTi(i,j),[1e-5 1000]);
    end
end

% pcolor(K,DTi,SIG); shading flat;
% contour(K,DTi,SIG); shading flat;

%% Computing local projections + high-end and low-end projections


% Mohsen has prepared the following mat file.
% load('/Users/ali/Downloads/T_kelvin_monthly_2000_2300.mat')
% load('/Users/ali/Downloads/T_kelvin_annual_2000_2100.mat')
% load('data/T_kelvin_annual_2000_2300_fixed.mat')
% T_proj_local = struct;
% T_proj_local.lat = Lat;
% T_proj_local.lon = Lon - 180;
% T_proj_local.TR = T_kelvin - repmat(T_kelvin(:,:,1), [1 1 size(T_kelvin, 3)]);
% 
% % T_proj_local.time_start = datenum(2000,6,30);
% % T_proj_local.time_end = datenum(2100,6,30);
% % T_proj_local.time = linspace(T_proj_local.time_start, T_proj_local.time_end, 101);
% 
% T_proj_local.time_start = datenum(2000,6,30);
% T_proj_local.time_end = datenum(2300,6,30);
% T_proj_local.time = linspace(T_proj_local.time_start, T_proj_local.time_end, 301);
% 
% save('Data/T_proj_local.mat', 'T_proj_local');
load('Data/T_proj_local.mat');

% load('Data/full-run.mat')
cd RCP
[ti, Tmax, Tmean, Tmin]=temperatureProjectionsCMIP5();
cd ..

[b,a] = butter(4,0.01,'low');
Tmin=filtfilt(b,a,Tmin);
Tmean=filtfilt(b,a,Tmean);
Tmax=filtfilt(b,a,Tmax);

Tmin(1:250) = interp1([1 250],[Tmin(1) Tmin(250)],1:1:250,'spline');
Tmax(1:250) = interp1([1 250],[Tmax(1) Tmax(250)],1:1:250,'spline');

Tmin(1:9) = Tmin(9)*ones(9,1);

ti_int = linspace(T_proj_local.time_start, T_proj_local.time_end, size(T_proj_local.TR, 3));
Tmax_ratio = 1+interp1(ti, (Tmax-Tmean)./Tmean, ti_int,'linear', 'extrap');
Tmin_ratio = 1+interp1(ti, (Tmin-Tmean)./Tmean, ti_int,'linear', 'extrap');
Tmin_ratio(Tmin_ratio<0.01) = 0.01;

Tmean_int = interp1(ti, Tmean, ti_int);

[X,Y] = meshgrid(T_proj_local.lat, T_proj_local.lon);

T_mean_stations = NaN(len,301);
T_low_stations = NaN(len,301);
T_high_stations = NaN(len,301);

for i=1:len
    id=ID(i);
    i
    % https://stackoverflow.com/questions/32481804/find-the-nearest-point-to-x-y-coordinates-on-mesh-matlab
    d = (S(id).lat-X).^2+(S(id).lon-Y).^2;
    [~, ind] = min(d(:));
    % resultX = X(ind) %// use that index to obtain the result
    % resultY = Y(ind)
    [ind_r, ind_c] = ind2sub(size(X), ind);
    % X(ind_r, ind_c)
    % Y(ind_r, ind_c)
    % The Temperature projection linked to the location:
    
    Ttmp = reshape(T_proj_local.TR(ind_r, ind_c, :), 1, []);
    Ttmp2 = Ttmp;

%     N=30;
%     Ttmp2=[flipud(Ttmp(1:N)), Ttmp, flipud(Ttmp(end-N:end))];
%     [b,a] = butter(2,0.02,'low');

    [b,a] = butter(4,0.01,'low'); %0.01
    Tsmooth=filtfilt(b,a,Ttmp2);
%     Tsmooth=Tsmooth(N+1:end-N-1);
   
    Tsmooth_shift = 0 - Tsmooth(1);
    Tsmooth = Tsmooth + Tsmooth_shift;

%     figure(13)
%     plot(Ttmp2); hold on
%     plot(Tsmooth)

%     plot(Tsmooth); hold on
%     drawnow
%     for ii = 1:301
%         if Tsmooth(ii)<0
%             pause(1)
%             break
%         end
%     end

    S(id).TR_raw = reshape(T_proj_local.TR(ind_r, ind_c, :), 1, []);
    S(id).TR_smooth = Tsmooth;
%     S(id).TR = S(id).TR_raw + Tsmooth_shift;
    S(id).TR = Tsmooth;
    assert(~any(isnan(S(id).TR)))
%     S(id).TR_high_end = S(id).TR.*(1+Tmax_ratio);
%     S(id).TR_low_end = S(id).TR.*(1+Tmin_ratio);
    S(id).TR_high_end = S(id).TR.*(Tmax_ratio);
    S(id).TR_low_end = S(id).TR.*(Tmin_ratio);
    
    S(id).TR_scenario = [S(ID(i)).TR_low_end', S(ID(i)).TR', S(ID(i)).TR_high_end'];
    S(id).ti = ti_int;


    T_mean_stations(i,:) = S(id).TR;
    T_low_stations(i,:) = S(id).TR.*(Tmin_ratio);
    T_high_stations(i,:) = S(id).TR.*(Tmax_ratio);

end

% ssss = 1050;
% figure(14)
% plot(T_mean_stations(ssss,:),'-k'); hold on
% plot(T_low_stations(ssss,:),'-b'); 
% plot(T_high_stations(ssss,:),'-r'); 

% Tmean_local = reshape(mean(mean(T_proj_local.TR, 1,'omitnan'),2,'omitnan'), 1, []);
% Tmean_local = interp1(ti_int, Tmean_local, ti);


% Global projections: mean-max-min
% figure;
% plot(ti, Tmean); hold on;
% plot(ti_int, Tmean_int.*(1+Tmax_ratio));
% plot(ti_int, Tmean_int.*(1+Tmin_ratio));
% datetick('x','keeplimits','keepticks');
% legend('Tmean','Tmax_comp. by ratio','Tmin comp. by ratio')
% % axis([datenum(2000,1,1) datenum(2100,1,1) 0 6]);
% axis([datenum(2000,1,1) datenum(2300,1,15) 0 16]);
% 
% % Local projection:
% figure; hold on;
% plot(T_proj_local.time, S(id).TR);
% plot(T_proj_local.time, S(id).TR_smooth);
% plot(T_proj_local.time, S(id).TR_high_end);
% plot(T_proj_local.time, S(id).TR_low_end);
% 
% datetick('x','keeplimits','keepticks');
% legend('Local TR', 'Local TR smooth', 'Local High-end TR', 'Local Low-end TR')


%% Calculations for ENSO vs. nonstationary GEV parameters

load('ENSO_ONI.mat');

ENSO_thresh = 0.833; %threshold to choose ENSO years based on ENSO index value (0.833 is the 0.5*std(ONI) value)
ENSO_yr = [];
ENSO_index = [];

for ff = 1:size(ENSO_ONI,1)
    if max(ENSO_ONI(ff,2:end)) > ENSO_thresh
        ENSO_yr = [ENSO_yr ENSO_ONI(ff,1)];
        ENSO_index = [ENSO_index max(ENSO_ONI(ff,2:end))];
    end
end

% ENSO_yr = [1952 1958 1964 1966 1969 1973 1983 1987 1988 1992 1995 1998 2003 2010 2016];
% ENSO_index = [1.2 1.8 1.4 2 1.1 2.1 2.2 1.2 1.7 1.7 1.1 2.4 1.3 1.6 2.6]; %ONI

NS_groups = [4 5 6 7 8 9]; %number of time periods for nonstationary analysis
start_yr = 1951;
end_yr = 2017;
max_events = 3; %top n maxima per year

for k = 1:length(NS_groups)
    NS_yr_lim = floor(linspace(start_yr,end_yr,NS_groups(k)+1));
    ENSO_index_avg = NaN(1,NS_groups(k));
    yearss = NaN(2,NS_groups(k));
    yearss_label = cell(1,NS_groups(k));
    for pp = 1:NS_groups(k)
        t_1 = datenum(append('01-Jun-',num2str(NS_yr_lim(pp))),'dd-mmm-yyyy');
        t_2 = datenum(append('01-Jun-',num2str(NS_yr_lim(pp+1))),'dd-mmm-yyyy');
        bin = [];
        for jj = 1:length(ENSO_yr)
            t_t = datenum(append('01-Jan-',num2str(ENSO_yr(jj))),'dd-mmm-yyyy');
            if t_t>t_1 && t_t<t_2
                bin = [bin ENSO_index(jj)];
            end
        end
        ENSO_index_avg(1,pp) = nanmean(bin);

        Yearss(1,pp) = NS_yr_lim(pp);
        Yearss(2,pp) = NS_yr_lim(pp+1);
        yearss_label(1,pp) = cellstr(append(num2str(nanmean(bin),'%.3f'),' ','(',num2str(NS_yr_lim(pp)),'-',num2str(NS_yr_lim(pp+1)),')'));

    end

    mu_st = NaN(1,len);
    sigma_st = NaN(1,len);
    k_st = NaN(1,len);
    mu_nonst = NaN(NS_groups(k),len);
    sigma_nonst = NaN(NS_groups(k),len);
    k_nonst = NaN(NS_groups(k),len);

    for i = 1:len

        [k i]
        
        max_temps = [];
        max_times = [];

        for j = 1:size(S(ID(i)).h_events,1)
            [bin_h,ind] = sort(S(ID(i)).h_events(j,:),'descend');
            bin_t = S(ID(i)).t_events(j,:);
            bin_t = bin_t(ind);
            max_temps = [max_temps bin_h(1:max_events)];
            max_times = [max_times bin_t(1:max_events)];
        end

        parmhat = gevfit(max_temps);
        mu_st(i) = parmhat(3);
        sigma_st(i) = parmhat(2);
        k_st(i) = parmhat(1);

        for ii = 1:NS_groups(k)
            t_1 = datenum(append('01-Jun-',num2str(NS_yr_lim(ii))),'dd-mmm-yyyy');
            t_2 = datenum(append('01-Jun-',num2str(NS_yr_lim(ii+1))),'dd-mmm-yyyy');

            ind2 = find(max_times>=t_1 & max_times<=t_2);
            max_temps_NS = max_temps(ind2);
            max_times_NS = max_times(ind2);

            if length(max_temps_NS) > floor(0.5*(NS_yr_lim(ii+1)-NS_yr_lim(ii)))

                parmhat = gevfit(max_temps_NS);
                mu_nonst(ii,i) = parmhat(3);
                sigma_nonst(ii,i) = parmhat(2);
                k_nonst(ii,i) = parmhat(1);

            end
        end
    end

    mu_dev = mu_nonst-repmat(mu_st,NS_groups(k),1);
    sigma_dev = sigma_nonst-repmat(sigma_st,NS_groups(k),1);
    k_dev = k_nonst-repmat(k_st,NS_groups(k),1);

    [ENSO_index_avg,ind4] = sort(ENSO_index_avg,'ascend');
    yearss_label = yearss_label(ind4);
    mu_dev = mu_dev(ind4,:);
    sigma_dev = sigma_dev(ind4,:);
    k_dev = k_dev(ind4,:);

    if k ==1
        x_starts = 0.12;
        y_starts = 0.78;
        x_len = 0.79;
        y_len = 0.12;

        ONI = [];
        for iii = 1:size(ENSO_ONI,1)
            ONI = [ONI ENSO_ONI(iii,2:end)];
        end

        figure(111);
        set(gcf,'PaperpositionMode','auto','Position',[100 50 1000 1000],'color','w');

        subplot('Position',[x_starts y_starts x_len y_len]); hold on; box on;
        plot(ENSO_ONI(1,1):1/12:ENSO_ONI(end,1)+1-1/12,ONI,'k');
        xx = ENSO_ONI(1,1):1/12:ENSO_ONI(end,1)+1-1/12;
        yy = ONI;
        xx = interp1(1:length(xx),xx,1:0.1:length(xx));
        yy = interp1(1:length(yy),yy,1:0.1:length(yy));
        idxp = yy>=0;
        area(xx(idxp),yy(idxp),'FaceColor',[0.8500 0.3250 0.0980]);
        idxn = yy<=0;
        area(xx(idxn),yy(idxn),'FaceColor',[0 0.4470 0.7410]);
        yline(std(ONI),'--','color',[0.4940 0.1840 0.5560],'linewidth',1);
        yline(-std(ONI),'--','color',[0.4940 0.1840 0.5560],'linewidth',1);
        area(xx,std(ONI)*ones(size(xx)),'FaceColor','w','edgeColor','none','FaceAlpha',0.7);
        area(xx,-std(ONI)*ones(size(xx)),'FaceColor','w','edgeColor','none','FaceAlpha',0.7);
        xlabel('Year'); xlim([1950 2017]);
        ylabel('ONI'); ylim([-2.75 2.75]);
        set(gca,'FontSize',12,'FontName','Constantia');
        legend('','El Niño event','La Niña event','\pm \sigma_{ONI}','','','','Location','southwest','FontSize',10,'Box','on','orientation','horizontal');
        title('The evolution of ENSO (1950-2017)','FontSize',16,'FontName','Constantia');

        annotation(gcf,'textbox',...
            [0.120000000000001 0.904221311475419 0.020999999999999 0.0220081967213137],...
            'VerticalAlignment','middle',...
            'String','A',...
            'Interpreter','none',...
            'HorizontalAlignment','center',...
            'FontWeight','bold',...
            'FontSize',18,...
            'FontName','Nomada Incise',...
            'FitBoxToText','off',...
            'BackgroundColor',[0.8 0.8 0.8]);

    end

    xlims = [1.44 1.63; 1.35 1.85; 1.35 2.06; 1.04 2.06; 1.04 2.06; 0.84 2.23];

    subplot(4,2,k+2); hold on; box on; grid on
    boxplot(k_dev',ENSO_index_avg,'Positions',ENSO_index_avg,'PlotStyle','compact');
    hh = findobj(gcf,'tag','Outliers');
    set(hh,'MarkerSize',0.000001);
    medlinreg = polyfit(ENSO_index_avg,nanmedian(k_dev'),1);
    medlinval = polyval(medlinreg,ENSO_index_avg);
    plot(ENSO_index_avg,medlinval,'--','color',[0.8500 0.3250 0.0980],'linewidth',2)
    % xlabel('$\overline{ONI}$','Interpreter','latex'); 
    xlim(xlims(k,:))
    xticks(ENSO_index_avg);
    xticklabels(string(yearss_label)); xtickangle(20)
    ylabel('$\Delta$k','Interpreter','latex'); 
    ylim([-0.69 0.69])
    set(gca,'FontSize',12,'FontName','Constantia');
    title(append('$\Delta$k vs. $\overline{ONI}$ [','y=',num2str(medlinreg(1),'%.3f'),'x+(',num2str(medlinreg(2),'%.3f'),')]'),'FontSize',15,'Interpreter','latex');

    labelsss = {'A','B','C','D','E','F','G'};
    x_sepp = 0.458;
    y_sepp = 0.223;

    annotation(gcf,'textbox',...
        [0.131000000000001+(mod(k+1,2))*x_sepp 0.69364754098361-floor((k-1)/2)*y_sepp 0.0209999999999992 0.0220081967213137],...
        'VerticalAlignment','middle',...
        'String',string(labelsss(k+1)),...
        'Interpreter','none',...
        'HorizontalAlignment','center',...
        'FontWeight','bold',...
        'FontSize',18,...
        'FontName','Nomada Incise',...
        'FitBoxToText','off',...
        'BackgroundColor',[0.8 0.8 0.8]);

    annotation(gcf,'textbox',...
        [0.171+(mod(k+1,2))*x_sepp 0.517418032786886-floor((k-1)/2)*y_sepp 0.282 0.0256147540983606],...
        'VerticalAlignment','middle',...
        'String','$\overline{ONI}$ (Period)',...
        'Interpreter','latex',...
        'HorizontalAlignment','center',...
        'FontSize',12,...
        'FitBoxToText','off',...
        'EdgeColor','none');

end

exportgraphics(gcf,'figures/ENSO_effect_on_k.png','Resolution',300)

%%
stop
%% clustering

Nclusters=4;

mu=[S(ID).mu]; sigma=[S(ID).sigma]; k=[S(ID).k];
MHHW=[S(ID).summer_avg];

RECLUSTER=0;
if RECLUSTER==1
    idx = kmeans([((mu-nanmean(mu))/nanstd(mu))' ((sigma-nanmean(sigma))/nanstd(sigma))'...
        ((k-nanmean(k))/nanstd(k))'],Nclusters,'MaxIter',1000);
    save('clusters.mat','idx');
else
    load(fullfile('data', 'clusters.mat'));
end

% color based on cluster 
   
cmap=[0   255 0;   ...
   255 0   0; ...
   0.9290*255 0.6940*255 0.1250*255;   ... 
   135 206 250;   ...
   ]/255;   

% cmap=[0.4660 0.6740 0.1880;   ...
%    0.6350 0.0780 0.1840; ...
%    0.9290 0.6940 0.1250;   ... 
%    0 0.4470 0.7410;   ...
%    ];    

FontName = 'Alverata';


% figure(10);
% scatter3([S(ID).mu],[S(ID).sigma],[S(ID).k],30*ones(length([S(ID).lon]),1),idx,'filled');
% xlabel('mu');
% ylabel('sigma');
% zlabel('k');
% colormap(cmap);

%cmap=colormap(jet(Nclusters));
   
han1=figure(1); set(han1,'PaperpositionMode','auto','Position',[200 100 1000 800],'color','w'); hold on; box on;
subplot('position',[0.056 0.49 0.542 0.451666666666667]);
hold on; 
geoshow('landareas.shp','FaceColor',[0.6 0.6 0.6],'EdgeColor',[0.6 0.6 0.6]); box on; axis equal
% axis equal; axis tight; box on; shading flat; hold on
% scatter([S(ID).longitude],[S(ID).latitude],30*ones(length([S(ID).longitude]),1),idx,'filled');
% scatter([S(ID).longitude],[S(ID).latitude],30*ones(length([S(ID).longitude]),1),[S(ID).mu],'filled');
%scatter([S(ID).lon],[S(ID).lat],30*ones(length([S(ID).lon]),1),[S(ID).mu],'filled');
scatter([S(ID).lon],[S(ID).lat],30*ones(length([S(ID).lon]),1),idx,'filled');
axis([-180 -60 0 80]);
colormap(cmap); 
set(gca,'FontSize',14,'FontName',FontName);
title('North America climate stations','FontSize',18);

% Borders: https://www.mathworks.com/matlabcentral/fileexchange/50390-borders
cd borders_v2/borders
borders('countries','k','nomap','linewidth',0.25)
cd ../..


id=find(strcmp({S(ID).id}','USW00094846')); % station of interest - Chicago
han1=plot(S(ID(id)).lon,S(ID(id)).lat,'ko'); set(han1,'LineWidth',2)
annotation(gcf,'textarrow',[0.443 0.467],[0.795 0.735],...
    'String',{'Chicago'},'FontName',FontName);
% 
% annotation(gcf,'textbox',[0.175 0.6 0.128 0.0524999999999999],...
%     'String',{'GEV model details','shown in panel B'},...
%     'FitBoxToText','off','FontSize',9,'linestyle','none','color',[0.6 0.6 0.6]);
% annotation(gcf,'arrow',[0.181 0.159],[0.635 0.615],'color',[0.6 0.6 0.6]);


% SUBPLOT B
subplot('position',[0.68 0.485 0.3 0.455]); hold on; box on;

axis([ri 100 20 50])
set(gca,'FontSize',14,'Xscale','log','FontName',FontName);

S(ID(id)).h_max=sort(reshape(S(ID(id)).h_events(:,1:3),[],1),'descend');
Nmax=length(S(ID(id)).h_max);
% max=S(ID(id)).yearsWithData*3;
% Plot results
P1=(0.001:0.001:1)'; E1=1-P1;
E_data=linspace(1/(Nmax+1),Nmax/(Nmax+1),Nmax)';
P_data=1-E_data;
ri=1/3;
%Tr_data=S(ID(id)).ri./E_data;
Tr_data=ri./E_data;
%ri=S.ri;

Tr1=(0.01:0.01:150)';
han1=plot(Tr_data,S(ID(id)).h_max,'o');set(han1,'MarkerFaceColor',[255 0 0]/255,'MarkerEdgeColor',[255 0 0]/255,'MarkerSize',5);
han1=plot(Tr1,x_Tr(Tr1,mu(id),sigma(id),k(id)),'b'); set(han1,'LineWidth',2);

GW=x_Tr(50,mu(id),sigma(id),k(id))-x_Tr(1,mu(id),sigma(id),k(id)); % ???? Is this how we get global warming?

han1=plot(Tr1,x_Tr(Tr1,mu(id)+GW,sigma(id),k(id)),'r');
set(han1,'LineWidth',2);

plot([0.1 100],[x_Tr(50,mu(id),sigma(id),k(id)) x_Tr(50,mu(id),sigma(id),k(id))],'k--');

%plot([1 1],[0 2],'k--');

annotation(gcf,'arrow',[0.738 0.738], [0.72 0.792],'LineWidth',2);
% annotation(gcf,'arrow',[0.737 0.737], [0.61118 0.71878],'LineWidth',2);
%annotation(gcf,'textbox',[0.701 0.62875 0.045 0.03875],'String',{'GW'},'FitBoxToText','off','FontWeight','bold','LineStyle','none');
annotation(gcf,'textbox',[0.702 0.726 0.045 0.0387],'String',{'TR'},'FitBoxToText','off','FontWeight','bold','LineStyle','none','FontSize',16,'FontName',FontName);
% 
annotation(gcf,'arrow',[0.734 0.734],[0.60155 0.67495],'LineWidth',2);
annotation(gcf,'textbox',[0.698 0.61377 0.045 0.0387],'String',{'TR'},'FitBoxToText','off','FontWeight','bold','LineStyle','none','FontSize',16,'FontName',FontName);

han1=plot([0.1 100],[MHHW(id)+GW MHHW(id)+GW],'r--'); set(han1,'LineWidth',2);

TR=[1 50];
for ii=1:length(TR)
    han1=plot(TR(ii),x_Tr(TR(ii),mu(id),sigma(id),k(id)),'kx'); set(han1,'MarkerSize',10,'LineWidth',2);
end

%han1 = text(.79,.81, 'ali'); set(han1, 'FontSize', 16)
han1=text(2.77223050497508, 35.1557987627863,[num2str(TR(1)),'-yr event = ',num2str(x_Tr(TR(1),mu(id),sigma(id),k(id)),3),' \circC']); set(han1,'FontSize',16,'FontName',FontName);
han1=text(1.221, 46.737 ,[num2str(TR(2)),'-yr event = ',num2str(x_Tr(TR(2),mu(id),sigma(id),k(id)),3), ' \circC']); set(han1,'FontSize',16,'FontName',FontName);
annotation(gcf,'arrow',[0.8832 0.938],[0.87142 0.80125]);
annotation(gcf,'arrow',[0.787 0.748],[0.715 0.71625]);

xlabel('Return period [years]');
ylabel('Extreme temperature [\circC]');
set(gca,'FontSize',14,'Xscale','log','FontName',FontName);
title('GEV model for Chicago, Illinois','FontSize',18)

han1=plot([0.1 100],[MHHW(id) MHHW(id)],'k--'); set(han1,'LineWidth',2);
% han1=text(1.3,MHHW(id)+0.04,['GW = ',num2str(MHHW(id),3),' m']); set(han1,'FontSize',16);
han1=text(1.7241, 28.4955,['MST = ',num2str(MHHW(id),3),' \circC']); set(han1,'FontSize',16,'FontName',FontName);

han1=legend('Temperature data','GEV model','GEV model + TR'); set(han1,'location','southeast','FontSize',16);

axis([ri 100 20 50])
%axis([ri 100 ylim])

width=0.24;
height=1.25*width;
xsep=0.09;
% SUBPLOT C
subplot('position',[0.085-0.01 0.11 width height]); hold on; box on;
for i=1:Nclusters
    han1=plot(mu(idx==i),sigma(idx==i),'o');set(han1,'MarkerFaceColor',cmap(i,:),'MarkerSize',8,'MarkerEdgeColor','k');
end
set(gca,'FontSize',14,'FontName',FontName)
xlabel('\mu','FontSize',18);
ylabel('\sigma','FontSize',18);
axis([10 50 0.5 3.5])


subplot('position',[0.085+(width+xsep)-0.01 0.11 width height]); hold on; box on;
for i=1:Nclusters
    han1=plot(mu(idx==i),k(idx==i),'o');set(han1,'MarkerFaceColor',cmap(i,:),'MarkerSize',8,'MarkerEdgeColor','k');
end
set(gca,'FontSize',14,'FontName',FontName)
xlabel('\mu','FontSize',18);
ylabel('$k$','Interpreter','latex','FontSize',18);
axis([10 50 -0.5 0.21])

subplot('position',[0.085+2*(width+xsep)-0.01 0.11 width height]); hold on; box on;
for i=1:Nclusters
    han1=plot(sigma(idx==i),k(idx==i),'o');set(han1,'MarkerFaceColor',cmap(i,:),'MarkerSize',8,'MarkerEdgeColor','k');
end
set(gca,'FontSize',14,'FontName',FontName)
xlabel('\sigma','FontSize',18);
ylabel('$k$','Interpreter','latex','FontSize',18);
axis([0.5 3.5 -0.5 0.21])
xticks([0.5:0.5:3.5])

annotation(gcf,'textbox',...
    [0.0630000000000006 0.49875 0.0349999999999997 0.04375],...
    'VerticalAlignment','middle',...
    'String','A',...
    'Interpreter','none',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',30,...
    'FontName','Nomada Incise',...
    'FitBoxToText','off',...
    'BackgroundColor',[0.8 0.8 0.8]);
annotation(gcf,'textbox',...
    [0.687000000000001 0.49875 0.0349999999999997 0.04375],...
    'VerticalAlignment','middle',...
    'String','B',...
    'Interpreter','none',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',30,...
    'FontName','Nomada Incise',...
    'FitBoxToText','off',...
    'BackgroundColor',[0.8 0.8 0.8]);
annotation(gcf,'textbox',...
    [0.0920000000000008-0.01 0.11875 0.0349999999999995 0.0437499999999998],...
    'VerticalAlignment','middle',...
    'String','C',...
    'Interpreter','none',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',30,...
    'FontName','Nomada Incise',...
    'FitBoxToText','off',...
    'BackgroundColor',[0.8 0.8 0.8]);
annotation(gcf,'textbox',...
    [0.422000000000001-0.01 0.11875 0.0349999999999995 0.0437499999999998],...
    'VerticalAlignment','middle',...
    'String','D',...
    'Interpreter','none',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',30,...
    'FontName','Nomada Incise',...
    'FitBoxToText','off',...
    'BackgroundColor',[0.8 0.8 0.8]);
annotation(gcf,'textbox',...
    [0.752000000000001-0.01 0.11875 0.0349999999999996 0.0437499999999998],...
    'VerticalAlignment','middle',...
    'String','E',...
    'Interpreter','none',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',30,...
    'FontName','Nomada Incise',...
    'FitBoxToText','off',...
    'BackgroundColor',[0.8 0.8 0.8]);


% print(gcf,'figures/Figure1_clusters.jpg','-djpeg','-r300');
exportgraphics(gcf,'figures/Figure2.png','Resolution',300)
% exportgraphics(gcf,'figures/test2.png','Resolution',300)


%% Plot Delta tempereture (Figure 2)

han1=figure(2); set(han1,'PaperpositionMode','auto','Position',[200 100 1000 500],'color','w');
subplot('position',[0.07 0.11 0.47 .82]); hold on; box on;
han1=plot([S(ID).k],DeltaT_50yr_1yr,'ko'); set(han1,'linewidth',2);
scatter([S(ID).k],DeltaT_50yr_1yr,40,[S(ID).sigma],'o','filled'); 
%caxis([0 4]); 
caxis([0.5 3.5]); 
colormap(flipud(jet)); 
han1=colorbar; set(han1,'visible','off')
% tmp=linspace(0.4,3.2,5); % BUG: has to be flipped
tmp=flip(linspace(0.4,3.2,5));
% contour(K,DTi,SIG,[6.4 3.2 1.6 0.8 0.4],'k--')

contour(K,DTi,SIG,tmp,'k--') % original
text(-0.1101, 10.213,['$\sigma =' num2str(tmp(1)) ' $'],'Interpreter','latex','Rotation',65,'FontSize',14,'color','k');
text(-0.0249, 10.154,['$\sigma =' num2str(tmp(2)) ' $'],'Interpreter','latex','Rotation',65,'FontSize',14,'color','k');
text(0.0817, 10.21,['$\sigma =' num2str(tmp(3)) ' $'],'Interpreter','latex','Rotation',65,'FontSize',14,'color','k');
text(0.1788, 8.4287,['$\sigma =' num2str(tmp(4)) ' $'],'Interpreter','latex','Rotation',62,'FontSize',14,'color','k');
text(0.1288, 2.7215,['$\sigma =' num2str(tmp(5)) ' $'],'Interpreter','latex','Rotation',35,'FontSize',14,'color','k');



LABELS=1;
if LABELS
    i = ID_cities;
    han1=plot([S(ID(i)).k],DeltaT_50yr_1yr(i),'ko'); set(han1,'linewidth',2);
    
%     annotation(gcf,'textarrow',[0.213 0.161],[0.138 0.196],...
%     'String',{'Honolulu'});
%     annotation(gcf,'textarrow',[0.233 0.266],[0.568 0.484],...
%     'String',{'San Francisco'});
%     annotation(gcf,'textarrow',[0.182 0.222],[0.5 0.408],...
%     'String',{'Seattle'});
%     annotation(gcf,'textarrow',[0.146 0.237],[0.426 0.364],...
%     'String',{'Juneau'});
%     annotation(gcf,'textarrow',[0.297 0.24],[0.164 0.31],...
%     'String',{'Dallas'});
%     annotation(gcf,'textarrow',[0.399 0.313],[0.408 0.398],...
%     'String',{'Chicago'});
%     annotation(gcf,'textarrow',[0.354 0.26],[0.186 0.312],...
%     'String',{'Boston'});
%     annotation(gcf,'textarrow',[0.358 0.271],[0.368 0.338],...
%     'String',{'Washington D.C.'});

    annotation(gcf,'textarrow',[0.213 0.161],[0.138 0.196],...
    'String',{'Honolulu'},'FontName',FontName);
    annotation(gcf,'textarrow',[0.233 0.275],[0.568 0.48],...
    'String',{'San Francisco'},'FontName',FontName);
    annotation(gcf,'textarrow',[0.182 0.23],[0.5 0.406],...
    'String',{'Seattle'},'FontName',FontName);
    annotation(gcf,'textarrow',[0.146 0.243],[0.426 0.366],...
    'String',{'Juneau'},'FontName',FontName);
    annotation(gcf,'textarrow',[0.297 0.249],[0.164 0.308],...
    'String',{'Dallas'},'FontName',FontName);
    annotation(gcf,'textarrow',[0.399 0.322],[0.408 0.398],...
    'String',{'Chicago'},'FontName',FontName);
    annotation(gcf,'textarrow',[0.354 0.268],[0.186 0.316],...
    'String',{'Boston'},'FontName',FontName);
    annotation(gcf,'textarrow',[0.358 0.28],[0.368 0.338],...
    'String',{'Washington D.C.'},'FontName',FontName);
    
end
axis([-0.55 .25 0 14]); % The next plot has to have the same axis
xlabel('$k$','FontSize',14,'interpreter','latex');
ylabel(['\DeltaT [\circC]'],'FontSize',14); 
set(gca,'FontSize',14,'XTick',[-1:0.25:1],'FontName',FontName);
title('50-yr event \rightarrow 1-yr event','FontSize',18);

% Right subplot 
subplot('position',[0.49 0.11 0.49 .82]); hold on; box on;
han1=plot([S(ID).k],DeltaT_1yr_SUMTEMP,'ko'); set(han1,'linewidth',2);
han2=scatter([S(ID).k],DeltaT_1yr_SUMTEMP,40,[S(ID).sigma],'o','filled'); 
caxis([0.5 3.5]); colormap(flipud(jet)); colorbar;

if LABELS
    i = ID_cities;
    
%     hax = [S(ID(i)).k]
%     hay = DeltaT_1yr_SUMTEMP(i)
    
    han1=plot([S(ID(i)).k],DeltaT_1yr_SUMTEMP(i),'ko'); set(han1,'linewidth',2)
    
%     annotation(gcf,'textarrow',[0.615990099009901 0.597990099009901],...
%     [0.186227180527383 0.240227180527383],'String',{'Honolulu'});
%     annotation(gcf,'textarrow',[0.807445544554456 0.721782178217822],...
%     [0.862170385395538 0.908722109533469],'String',{'San Francisco'});
%     annotation(gcf,'textarrow',[0.603960396039604 0.664009900990099],...
%     [0.860040567951318 0.721079107505071],'String',{'Seattle'});
%     annotation(gcf,'textarrow',[0.804980198019802 0.694980198019802],...
%     [0.746908722109533 0.654908722109533],'String',{'Juneau'});
%     annotation(gcf,'textarrow',[0.586 0.672],...
%     [0.384511156186613 0.428511156186613],'String',{'Dallas'});
%     annotation(gcf,'textarrow',[0.824990099009901 0.757990099009901],...
%     [0.604766734279919 0.560766734279919],'String',{'Chicago'});
%     annotation(gcf,'textarrow',[0.577990099009901 0.692990099009901],...
%     [0.556908722109533 0.630908722109533],'String',{'Boston'});
%     annotation(gcf,'textarrow',[0.8 0.715],...
%     [0.241379310344828 0.478653144016227],'String',{'Washington D.C.'});


    annotation(gcf,'textarrow',[0.615990099009901-0.015 0.6-0.015],...
    [0.186227180527383 0.24],'String',{'Honolulu'},'FontName',FontName);
    annotation(gcf,'textarrow',[0.807445544554456-0.015 0.731-0.015],...
    [0.862170385395538 0.906],'String',{'San Francisco'},'FontName',FontName);
    annotation(gcf,'textarrow',[0.603960396039604-0.015 0.673-0.015],...
    [0.860040567951318 0.718],'String',{'Seattle'},'FontName',FontName);
    annotation(gcf,'textarrow',[0.555 0.674],[0.626 0.65],'String',{'Juneau'},'FontName',FontName);
    annotation(gcf,'textarrow',[0.586-0.015 0.679-0.015],...
    [0.384511156186613 0.43],'String',{'Dallas'},'FontName',FontName);
    annotation(gcf,'textarrow',[0.824990099009901-0.015 0.766-0.015],...
    [0.604766734279919 0.56],'String',{'Chicago'},'FontName',FontName);
    annotation(gcf,'textarrow',[0.802 0.698],[0.718 0.65],'String',{'Boston'},'FontName',FontName);
    annotation(gcf,'textarrow',[0.8-0.019 0.724-0.019],...
    [0.241379310344828 0.48],'String',{'Washington D.C.'},'FontName',FontName);

end


axis([-0.55 .25 0 14]); % BUG: This has to be the same the other subplot
% axis([-0.55 .25 0 12]);
xlabel('$k$','FontSize',14,'interpreter','latex'); set(gca,'FontSize',14,'XTick',[-1:0.25:1],'FontName',FontName);
set(gca,'YTickLabel',[]);

title('1-yr event \rightarrow MST','FontSize',18,'FontName',FontName);

annotation(gcf,'textbox',[0.924 0.924 0.03 0.078],'String',{'$\sigma$'},...
    'LineStyle','none','Interpreter','latex','FontSize',20,'FitBoxToText','off');
annotation(gcf,'textbox',[0.935000000000001 0.922 0.0549999999999993 0.078],'String',{' [\circ C]'},...
    'LineStyle','none','FontSize',16,'FitBoxToText','off','FontName',FontName);


annotation(gcf,'textbox',...
    [0.0780000000000012 0.842750000000002 0.0349999999999995 0.07125],...
    'VerticalAlignment','middle',...
    'String','A',...
    'Interpreter','none',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',30,...
    'FontName','Nomada Incise',...
    'FitBoxToText','off',...
    'BackgroundColor',[0.8 0.8 0.8]);
annotation(gcf,'textbox',...
    [0.498000000000001 0.842750000000003 0.0349999999999996 0.07125],...
    'VerticalAlignment','middle',...
    'String','B',...
    'Interpreter','none',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',30,...
    'FontName','Nomada Incise',...
    'FitBoxToText','off',...
    'BackgroundColor',[0.8 0.8 0.8]);
% print(gcf,'figures/Figure2_delta_temperature.jpg','-djpeg','-r300');
exportgraphics(gcf,'figures/Figure3.png','Resolution',300)
% exportgraphics(gcf,'figures/test2.png','Resolution',300)


% -------------------------------------------------------------------------
% Statistics: 
% hist(DeltaT_50yr_1yr)
ftemp = @(x) sum(DeltaT_50yr_1yr<x & DeltaT_50yr_1yr>0)/length(DeltaT_50yr_1yr)*100
f_time_to_get_temp = @(x) datetime(ti(find( Tmean > x, 1 )),'ConvertFrom','datenum')
%plot(1:.1:10,ftemp(1:0.1:10))
disp(['50yr becoming y1' string(ftemp(5)), '% of stations'])
disp(['50yr becoming y1' string(ftemp(6.5)), '% of stations'])
disp(['DT of two plots - average:'  string(mean(DeltaT_1yr_SUMTEMP))])
disp(['DT of two plots - standard deviation:'  string(std(DeltaT_1yr_SUMTEMP))])
sum(DeltaT_1yr_SUMTEMP<5)/length(DeltaT_50yr_1yr)*100

% histogram(DeltaT_1yr_SUMTEMP-DeltaT_50yr_1yr)
%% plot temperature rise horizons plot - Figure 3



TIME1=NaN(len,1);
TIME2=NaN(len,1);

TIME1_max=NaN(len,1);
TIME1_min=NaN(len,1);

TIME2_max=NaN(len,1);
TIME2_min=NaN(len,1);


% RCP={'RCP 2.6','RCP 4.5','RCP 6.0','RCP 8.5'};   
s=4;
% s=1; % RCP 8.5
% s=2; % RCP 4.5
% s=3; % RCP 2.6

% TEMP=0:20;

% REMOVE=15;
% RCPt(end-REMOVE:end)=[];
% RCPTEMP(end-REMOVE:end,:)=[];
%TIMEEX=interp1(RCPTEMP(:,s), RCPt,TEMP,'linear','extrap');
% TIMEEX=interp1(Tmean, ti,TEMP,'linear','extrap'); % NEW Commented

% if 1
%     figure;
%     plot(ti,Tmax,'-b.')%,TIME1,TEMP,'ro');
%     stop
% end

for i=1:len 
    
    % find the times when SLR is the same as the Delta water level
%     TIME1(i)=interp1(RCPTEMP(:,s), RCPt,DeltaT_50yr_1yr(i),'linear','extrap');
    
%     TIME1(i)=interp1(Tmean, ti, DeltaT_50yr_1yr(i),'linear','extrap');
%     TIME1_max(i)=interp1(Tmin, ti, DeltaT_50yr_1yr(i),'linear','extrap');
%     TIME1_min(i)=interp1(Tmax, ti, DeltaT_50yr_1yr(i),'linear','extrap');
%   NEW changes: using local projections:
    TIME1(i)=interp1(S(ID(i)).TR, T_proj_local.time, DeltaT_50yr_1yr(i),'linear','extrap');
    TIME1_max(i)=interp1(S(ID(i)).TR_low_end, T_proj_local.time, DeltaT_50yr_1yr(i),'linear','extrap');
    TIME1_min(i)=interp1(S(ID(i)).TR_high_end, T_proj_local.time, DeltaT_50yr_1yr(i),'linear','extrap');
    
    
%     TIME2(i)=interp1(Tmean, ti,DeltaT_1yr_SUMTEMP(i),'linear','extrap');
% NEW change: using local projections
    TIME2(i)=interp1(S(ID(i)).TR, T_proj_local.time,DeltaT_1yr_SUMTEMP(i),'linear','extrap');
    
%     figure
%     plot(TIMEEX, TEMP,'r-',RCPt,RCPTEMP(:,s),'b-',TIME1(i),DeltaT_50yr_1yr(i),'ro',[datenum(2000,1,1) datenum(2150,1,1)],[DeltaT_50yr_1yr(i) DeltaT_50yr_1yr(i)],'k--',[TIME1(i) TIME1(i)],[-1 7],'k--'); datetick('x','keeplimits');
%     
%     stop
%     
end

% Plot time horizon

YLIM1=datenum(2000,1,1);
YLIM2=datenum(2300,1,1);

han1=figure(3); set(han1,'PaperpositionMode','auto','Position',[200 100 1000 500],'color','w');
subplot('position',[0.07 0.11 0.41 .82]); hold on; box on;
%contour(K,TIME,SIG,[0.02 0.05 0.1 0.2],'k--');
han1=scatter([S(ID).k],TIME1,40,[S(ID).sigma],'o','filled'); 
%caxis([0 4]); colormap(jet); %colorbar;
caxis([0.5 3.5]); colormap(flipud(jet)); %colorbar;

yneg = abs(TIME1-TIME1_min);
ypos = abs(TIME1-TIME1_max);

xneg = zeros(size(TIME1));
xpos = zeros(size(TIME1));

ERRORBARS=0;
if ERRORBARS==1
    han1=errorbar([S(ID).k],TIME1,yneg,ypos,xneg,xpos,'.k'); set(han1,'color',[0.666 0.666 0.666])
end

han1=plot([S(ID).k],TIME1,'ko'); set(han1,'linewidth',2);
han1=scatter([S(ID).k],TIME1,40,[S(ID).sigma],'o','filled'); 
% caxis([0 0.2]); colormap(jet); %colorbar;
axis([-0.55 .25 YLIM1 YLIM2]);
LABELS=1;
if LABELS
    
    han1=plot([S(ID(ID_cities)).k],TIME1(ID_cities),'ko'); set(han1,'linewidth',2);
    
    % ADD labels
%     a=axis;
%     p=han1.Parent.Position;
%     x=@(x) (x-a(1))/(a(2)-a(1))-p(1);
%     y=@(y) (y-a(3))/(a(4)-a(3))-p(2);

    annotation(gcf,'textarrow',[0.197 0.163],[0.16 0.194],'String',cities(1,2),'FontName',FontName);
    annotation(gcf,'textarrow',[0.242 0.275],[0.538 0.414],'String',cities(2,2),'FontName',FontName);
    annotation(gcf,'textarrow',[0.148 0.229],[0.42 0.314],'String',cities(3,2),'FontName',FontName);
    annotation(gcf,'textarrow',[0.181 0.246],[0.482 0.352],'String',cities(4,2),'FontName',FontName);
    annotation(gcf,'textarrow',[0.279 0.247],[0.172 0.282],'String',cities(5,2),'FontName',FontName);
    annotation(gcf,'textarrow',[0.415 0.32],[0.442 0.354],'String',cities(6,2),'FontName',FontName);
    annotation(gcf,'textarrow',[0.33 0.266],[0.18 0.306],'String',cities(7,2),'FontName',FontName);
    annotation(gcf,'textarrow',[0.364 0.279],[0.272 0.34],'String',cities(8,2),'FontName',FontName);

end

set(gca,'FontSize',14,'XTick',[-1:0.25:1],'YTick',datenum(2000:25:2300,1,1),'FontName',FontName); 
datetick('y','keeplimits','keepticks');

xlabel('$k$','FontSize',14,'interpreter','latex','interpreter','latex'); ylabel('Year','FontSize',14);
title('50-yr event \rightarrow 1-yr event','FontSize',18);

subplot('position',[0.49 0.11 0.5 .82]); hold on; box on;
han1=scatter([S(ID).k],TIME2,40,[S(ID).sigma],'o','filled'); 
% caxis([0 0.2]); colormap(jet); colorbar;

% yneg = abs(TIME2-TIME2_min);
% ypos = abs(TIME2-TIME2_max);
% 
% xneg = zeros(size(TIME2));
% xpos = zeros(size(TIME2));
% 
% ERRORBARS=1;
% if ERRORBARS;
%     han1=errorbar([S(ID).k],TIME2,yneg,ypos,xneg,xpos,'.k');  set(han1,'color',[0.666 0.666 0.666])
% end

han1=plot([S(ID).k],TIME2,'ko'); set(han1,'linewidth',2);
han1=scatter([S(ID).k],TIME2,40,[S(ID).sigma],'o','filled'); 
%caxis([0 4]); colormap(jet); colorbar;
caxis([0.5 3.5]); colormap(flipud(jet)); colorbar;
% axis([-0.55 .25 YLIM1 YLIM2]);
xlabel('$k$','FontSize',14,'interpreter','latex'); 
set(gca,'FontSize',14,'XTick',[-1:0.25:1],'YTick',datenum(2000:25:2300,1,1),'FontName',FontName);
if LABELS
    
    i = ID_cities;
    han1=plot([S(ID(i)).k],TIME2(i),'ko'); set(han1,'linewidth',2);
    
    % ADD labels
%     a=axis;
%     p=han1.Parent.Position;
%     x=@(x) (x-a(1))/(a(2)-a(1))-p(1);
%     y=@(y) (y-a(3))/(a(4)-a(3))-p(2);

%     annotation(gcf,'textarrow',[0.65 0.632],[0.228 0.282],'String',cities(1,2));
% %     annotation(gcf,'textarrow',[0.276 0.301],[0.616 0.536],'String',cities(2,2));
% %     annotation(gcf,'textarrow',[0.22 0.252],[0.486 0.44],'String',cities(3,2));
%     annotation(gcf,'textarrow',[0.633 0.704],[0.842 0.886],'String',cities(4,2));
%     annotation(gcf,'textarrow',[0.621 0.698],[0.398 0.476],'String',cities(5,2));
%     annotation(gcf,'textarrow',[0.84 0.773],[0.69 0.646],'String',cities(6,2));
%     annotation(gcf,'textarrow',[0.796 0.727],[0.882 0.862],'String',cities(7,2));
%     annotation(gcf,'textarrow',[0.817 0.734],[0.244 0.536],'String',cities(8,2));
    
    
    annotation(gcf,'textarrow',[0.605 0.587],[0.174 0.228],'String',cities(1,2),'FontName',FontName);
%     annotation(gcf,'textarrow',[0.276 0.301],[0.616 0.536],'String',cities(2,2),'FontName',FontName);
    annotation(gcf,'textarrow',[0.592 0.661],[0.658 0.53],'String',cities(3,2),'FontName',FontName);
    annotation(gcf,'textarrow',[0.639 0.681],[0.768 0.54],'String',cities(4,2),'FontName',FontName);
    annotation(gcf,'textarrow',[0.582 0.669],[0.334 0.38],'String',cities(5,2),'FontName',FontName);
    annotation(gcf,'textarrow',[0.84 0.757000000000001],[0.48 0.44],'String',cities(6,2),'FontName',FontName);
    annotation(gcf,'textarrow',[0.791 0.701],[0.606 0.546],'String',cities(7,2),'FontName',FontName);
    annotation(gcf,'textarrow',[0.725 0.709],[0.22 0.434],'String',cities(8,2),'FontName',FontName);
end

annotation(gcf,'textbox',[0.931 0.926 0.0300000000000001 0.0780000000000002],'String',{'$\sigma$'},...
    'LineStyle','none','Interpreter','latex','FontSize',20,'FitBoxToText','off');
annotation(gcf,'textbox',[0.948000000000001 0.924 0.0459999999999989 0.0780000000000002],'String',{'[\circ C]'},...
    'LineStyle','none','FontSize',16,'FitBoxToText','off','FontName',FontName);

axis([-0.55 .25 YLIM1 YLIM2]);
xlabel('$k$','FontSize',14,'interpreter','latex'); 
set(gca,'FontSize',14,'XTick',[-1:0.25:1],'YTick',datenum(2000:25:2300,1,1),'FontName',FontName); 
datetick('y','keeplimits','keepticks');
set(gca,'YTickLabel',[]);

title('1-yr event \rightarrow MST','FontSize',18,'FontName',FontName);

annotation(gcf,'textbox',...
    [0.0780000000000012 0.842750000000002 0.0349999999999995 0.07125],...
    'VerticalAlignment','middle',...
    'String','A',...
    'Interpreter','none',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',30,...
    'FontName','Nomada Incise',...
    'FitBoxToText','off',...
    'BackgroundColor',[0.8 0.8 0.8]);
annotation(gcf,'textbox',...
    [0.498000000000001 0.842750000000003 0.0349999999999996 0.07125],...
    'VerticalAlignment','middle',...
    'String','B',...
    'Interpreter','none',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',30,...
    'FontName','Nomada Incise',...
    'FitBoxToText','off',...
    'BackgroundColor',[0.8 0.8 0.8]);

% print(gcf,'figures/Figure3_time_horizons.jpg','-djpeg','-r300');
exportgraphics(gcf,'figures/Figure4.png','Resolution',300)

%% Exceedance Plot with temperature (Appendix 1)

han1=figure(4); set(han1,'PaperpositionMode','auto','Position',[200 100 1200 600],'color','w');

GW=(0:0.1:10.0)';
    
E_SLR=NaN(length(GW),len);
TR_SLR=NaN(length(GW),len);

E_SLR_DATA=NaN(length(GW),len);
TR_SLR_DATA=NaN(length(GW),len);

E_SLR_GEV=NaN(length(GW),len);
TR_SLR_GEV=NaN(length(GW),len);
    
TR=50;

for i=1:1:len
% for i=1:30:len %  not all the stations are being counted 
    
    id=ID(i); % get the specific station
    
%     E0=S(id).ri/TR; % get the exceedance probability of the TR-year event
    E0=ri/TR;
    % GEV
    mu=S(id).mu;
    sig=S(id).sigma;
    k=S(id).k;

    % empirical exceedace probability of all data
    [Etmp,ia,~] = unique(S(id).Exceedance,'last');
    xtmp=S(id).x(ia);

    x0=x_Tr(TR,mu,sig,k); % find the Temperature corresponding to the 50-year event
    

    
    Ei=interp1(xtmp,Etmp,x0,'linear');
    %Ei=interp1(xtmp,Etmp,x0,'pchip');
    

    % empirical exceedace probability of the data of extreme events
    Nmax=3;
    t_max=S(id).t_events(:,1:Nmax);t_max=t_max(:);
    h_max=S(id).h_events(:,1:Nmax);h_max=h_max(:);
    [h_max,id]=sort(h_max,'descend'); 
    t_max=t_max(id);
    
%     len_data=length(S(id).h_max);
    len_data=length(h_max);
    E_data=linspace(1/(len_data+1),len_data/(len_data+1),len_data)';
    x_data=h_max;
    
    x0_data=interp1(E_data,x_data,E0,'linear','extrap');
    [~,ind]=unique(x_data); % ??? I added
    Ei_data=interp1(x_data(ind),E_data(ind),x0,'linear','extrap'); % ??? Unique!!
    
    x1=(0:0.001:2*max(x_data))';

    for kk=1:length(GW)
        
        % find empirical exceedance probability increase of the all values
        E_SLR(kk,i)=interp1(xtmp+GW(kk),Etmp,x0);
        TR_SLR(kk,i)=ri/E_SLR(kk);
        
        % find empirical exceedance probability increase of the extremes
        tmp=x_data+GW(kk);
        [~,ind]=unique(tmp); % ??? I added
        
        %E_SLR_DATA(kk,i)=interp1(tmp(ind),E_data(ind),x0_data,'linear','extrap');
        
        % find empirical exceedance probability increase of the extremes
        if GW(kk)>2
            E_SLR_DATA(kk,i)=interp1(tmp(ind),E_data(ind),x0_data);
        else
            E_SLR_DATA(kk,i)=interp1(tmp(ind),E_data(ind),x0_data,'linear','extrap');
        end
        
        TR_SLR_DATA(kk,i)=ri/E_SLR_DATA(kk,i);
        
        % find GEV probability increase
        E_SLR_GEV(kk,i)=E(x0,mu+GW(kk),sig,k);
        TR_SLR_GEV(kk,i)=ri/E_SLR_GEV(kk,i);
        
        PLOT=0;
        if PLOT
            
            figure;
            subplot(1,3,1);

            plot(xtmp,Etmp,'-b.',xtmp+GW(kk),Etmp,'-r.',...
                [0 max(x_data)+GW(end)],[Ei Ei],'k--',...
                x0,Ei,'mo',[x0 x0],[0 1],'k--',...
                x0,E_SLR(kk,i),'go',...
                [0 max(x_data)+GW(end)],[E_SLR(kk,i) E_SLR(kk,i)],'k--');
            
            axis([0 max(x_data)+GW(end) 0 1]);

            xlabel('Temp [C]')
            ylabel('Exceedance Prob.')

            subplot(1,3,2);
            
            plot(x_data,E_data,'-b.',x_data+GW(kk),E_data,'-r.',...
                [0.9*min(x_data) max(x_data)+GW(end)],[ri/TR ri/TR],'k--',...
                x0_data,E0,'mo',[x0_data x0_data],[0 1],'k--',...
                x0_data,E_SLR_DATA(kk,i),'go',...
                [0.9*min(x_data) max(x_data)+GW(end)],[E_SLR_DATA(kk,i) E_SLR_DATA(kk,i)],'k--',[x0 x0],[0 1],'g--');
            
            axis([0.9*min(x_data) max(x_data)+GW(end) 0 1]);
      
            subplot(1,3,3);
            
            plot(x1,E(x1,mu,sig,k),'-b.',x1,E(x1,mu+GW(kk),sig,k),'-r.',...
                [0.9*min(x_data) max(x_data)+GW(end)],[ri/TR ri/TR],'k--',...
                x_Tr(TR,mu,sig,k),E0,'mo',[x_Tr(TR,mu,sig,k) x_Tr(TR,mu,sig,k)],[0 1],'k--',...
                x_Tr(TR,mu,sig,k),E_SLR_GEV(kk,i),'go',...
                [0.9*min(x_data) max(x_data)+GW(end)],[E_SLR_GEV(kk,i) E_SLR_GEV(kk,i)],'k--');
            
            axis([0.9*min(x_data) max(x_data)+GW(end) 0 1]);
           
            xlabel('Temp [C]')
            ylabel('Exceedance Prob.')
            drawnow;
            stop

        end
        
    end

    color1=cmap(idx(i),:); % color based on clusters
%    color1=cmap(1,:); % color based on clusters
    
    %color1='b';           % single color
    
    ALPHA=0.15;
    
    xp0=0.09;
    yp0=0.12;
    width=0.28;
    height=0.77;
    xsep=0.02;
    
    subplot('position',[xp0 yp0 width height]); hold on; box on;
    han1=plot(GW,E_SLR(:,i)/Ei,'Color',color1); han1.Color(4)=ALPHA;  % set transparency
    set(han1,'LineWidth',2);
    title({'All data'},'Fontsize',14,'FontName',FontName);
    
    subplot('position',[xp0+(width+xsep) yp0 width height]); hold on; box on;
    han1=plot(GW,E_SLR_DATA(:,i)/E0,'Color',color1); han1.Color(4)=ALPHA;  % set transparency
    set(han1,'LineWidth',2);
    title({'Extreme events data'},'Fontsize',14,'FontName',FontName);
    
    subplot('position',[xp0+2*(width+xsep) yp0 width height]); hold on; box on;
    han1=plot(GW,E_SLR_GEV(:,i)/E0,'Color',color1); han1.Color(4)=ALPHA;  % set transparency
    set(han1,'LineWidth',2);
    title({'GEV model from extreme events data'},'Fontsize',14,'FontName',FontName);
end

T_max=max(GW);

for i=1:3
    subplot('position',[xp0+(i-1)*(width+xsep) yp0 width height]);
    plot(GW,2.^(GW/1),'--','LineWidth',3,'Color','k'); hold on
    plot(GW,2.^(GW/0.25),'--','LineWidth',3,'Color',[0.4940 0.1840 0.5560])
    plot(GW,2.^(GW/0.5),'--','LineWidth',3,'Color','b')
    plot(GW,2.^(GW/1.5),'--','LineWidth',3,'Color','r')
    
    plot([0 T_max],[150 150],'--k');
    
    set(gca,'FontSize',14);
    xlabel('\mu_T [\circC]','FontSize',16,'FontName',FontName);
    if i==1
        ylabel('$\frac{E}{E_0}$','FontSize',26,'Interpreter','Latex'); 
        set(gca,'Yscale','log','ytick',[1 2 4 8 16 32 64 128 256],'YMinorTick','off','Xtick',0:1:T_max,'FontName',FontName);
    else
        set(gca,'Yscale','log','ytick',[1 2 4 8 16 32 64 128 256],'yticklabel',[],'YMinorTick','off','Xtick',0:1:T_max,'FontName',FontName);
    end
    axis([0 T_max 1 256]);
    
    text(1.00684527839933, 61.5049107158222,'0.25 \circC','Rotation',80,'FontSize',14,'color',[0.4940 0.1840 0.5560],'FontName',FontName);
    text(2.48303569498516, 57.231344506145,'0.5 \circC','Rotation',75,'FontSize',14,'color','b','FontName',FontName);
    text(4.80654761727367, 41.78037710968,'1 \circC','Rotation',62,'FontSize',14,'color','k','FontName',FontName);
    text(6.59196436405182, 28.0478757570636,'1.5 \circC','Rotation',50,'FontSize',14,'color','r','FontName',FontName);

end

annotation(gcf,'textbox',[0.250000000000002 0.920000000000007 0.997500000000009 0.0633333333333341],'String',{'Normalized increase in the probability of exceedance (i.e., "doubling")'},'FitBoxToText','off','linestyle','none','FontWeight','bold','FontSize',18,'FontName',FontName);

annotation(gcf,'textbox',...
    [0.333833333333334 0.135000000000001 0.0286666666666665 0.0633333333333344],...
    'VerticalAlignment','middle',...
    'String',{'A'},...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',30,...
    'FontName','Nomada Incise',...
    'FitBoxToText','off',...
    'BackgroundColor',[0.8 0.8 0.8]);
annotation(gcf,'textbox',...
    [0.633833333333341 0.135000000000003 0.0286666666666665 0.0633333333333345],...
    'VerticalAlignment','middle',...
    'String','B',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',30,...
    'FontName','Nomada Incise',...
    'FitBoxToText','off',...
    'BackgroundColor',[0.8 0.8 0.8]);
annotation(gcf,'textbox',...
    [0.933833333333339 0.135000000000002 0.0286666666666664 0.0633333333333345],...
    'VerticalAlignment','middle',...
    'String','C',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',30,...
    'FontName','Nomada Incise',...
    'FitBoxToText','off',...
    'BackgroundColor',[0.8 0.8 0.8]);

% print(gcf,'figures/Figure4_rate_of_increase_with_GTR_Appendix1.jpg','-djpeg','-r300');
exportgraphics(gcf,'figures/FigureA2.png','Resolution',300)

%% Exceedance Plot with temperature 2 (Figure 4) (Appendix 2)

han1=figure(444); set(han1,'PaperpositionMode','auto','Position',[200 100 650 600],'color','w'); hold on; box on;

han1=figure(44); set(han1,'PaperpositionMode','auto','Position',[200 100 1200 600],'color','w');

GW=(0:0.1:10.0)'; 
% We don't have a global warming senario here because it's independant of
% time.
% GW=(0:0.01:1.0)'; % TODO: ???? What does GW should be?

% GW_threshold=0.25;        % SLR threshold for fit of sigma 2x
GW_threshold=3; % NEW (Sean rec.)
idsel=GW<GW_threshold;    % TODO: ???? What does this threshold should be?

EI=NaN(1,len);
    
E_SLR=NaN(length(GW),len);
TR_SLR=NaN(length(GW),len);
ODDS_SLR=NaN(length(GW),len);

E_SLR_DATA=NaN(length(GW),len);
TR_SLR_DATA=NaN(length(GW),len);
ODDS_SLR_DATA=NaN(length(GW),len);

E_SLR_GEV=NaN(length(GW),len);
TR_SLR_GEV=NaN(length(GW),len);
ODDS_SLR_GEV=NaN(length(GW),len);   

lat=NaN(len,1);
lon=NaN(len,1);

clus=NaN(len,1);

sigma_2x=NaN(len,1);
sigma_2x_DATA=NaN(len,1);
sigma_2x_GEV=NaN(len,1);

sigma_2x_fit=NaN(len,1);
sigma_2x_fit_DATA=NaN(len,1);
sigma_2x_fit_GEV=NaN(len,1);

TR=50;

TR_max=max(GW);

nnnn = 0; %no. of lines going out of the light dashed lines
% for i=1:1:len 
for i=1:1:len % I'm not covering all the stations
    i
    
    id=ID(i); % get the specific station
    
    E0=S(id).ri/TR; % get the exceedance probability of the TR-year event
%     E0=ri/TR; % The same as above since ri is constant. ???? Why ri is
%     constant

    % GEV
    mu=S(id).mu;
    sig=S(id).sigma;
    k=S(id).k;
    
    lat(i)=S(id).lat;
    lon(i)=S(id).lon;
    
    % cluster
    clus(i)=idx(i);

    % empirical exceedace probability of all data
    [Etmp,ia,~] = unique(S(id).Exceedance,'last');
    xtmp=S(id).x(ia);

    x0=x_Tr(TR,mu,sig,k); % find the Temperature corresponding to the 50-year event
    
    Ei=interp1(xtmp,Etmp,x0,'linear');
    %Ei=interp1(xtmp,Etmp,x0,'pchip');
    EI(i)=Ei;

%     % empirical exceedace probability of the data of extreme events
%     Nmax=3;
%     t_max=S(id).t_events(:,1:Nmax);t_max=t_max(:);
%     h_max=S(id).h_events(:,1:Nmax);h_max=h_max(:);
%     [h_max,id]=sort(h_max,'descend'); 
%     t_max=t_max(id);
%     
% %     len_data=length(S(id).h_max);
%     len_data=length(h_max);
%     E_data=linspace(1/(len_data+1),len_data/(len_data+1),len_data)';
%     x_data=h_max;
%     
%     x0_data=interp1(E_data,x_data,E0,'linear','extrap');
%     [~,ind]=unique(x_data); % ??? I added
%     Ei_data=interp1(x_data(ind),E_data(ind),x0,'linear','extrap'); % ??? Unique!!
%     
%     x1=(0:0.001:2*max(x_data))';

    % empirical exceedace probability of all data
    [Etmp,ia,~] = unique(S(id).Exceedance,'last');
    xtmp=S(id).x(ia);

    x0=x_Tr(TR,mu,sig,k); % find the water level corresponding to the 50-year event
    
    Ei=interp1(xtmp,Etmp,x0,'linear');
    EI(i)=Ei;
    
    % empirical exceedace probability of the data of extreme events
    len_data=length(S(id).h_max);
    E_data=linspace(1/(len_data+1),len_data/(len_data+1),len_data)';
    x_data=S(id).h_max;

    x0_data=interp1(E_data,x_data,E0,'linear','extrap');
    %Ei_data=interp1(x_data,E_data,x0,'linear');
    
    x1=(0:0.001:2*max(x_data))';
    

    for kk=1:length(GW)
        
        % find empirical exceedance probability increase of the all values
        E_SLR(kk,i)=interp1(xtmp+GW(kk),Etmp,x0);
        TR_SLR(kk,i)=ri/E_SLR(kk);
        
        % find empirical exceedance probability increase of the extremes
        [~,ind]=unique(x_data); % TODO: ???? I added. Is this Okay?
        
        % find empirical exceedance probability increase of the extremes
        if GW(kk)>2 % TODO: ???? Should this if statement be 2?
            E_SLR_DATA(kk,i)=interp1(x_data(ind) + GW(kk), E_data(ind),x0_data);
        else
            E_SLR_DATA(kk,i)=interp1(x_data(ind) + GW(kk), E_data(ind),x0_data,'linear','extrap');
        end
        TR_SLR_DATA(kk,i)=ri/E_SLR_DATA(kk,i);
        
        % find GEV probability increase
        E_SLR_GEV(kk,i)=E(x0,mu+GW(kk),sig,k);
        TR_SLR_GEV(kk,i)=ri/E_SLR_GEV(kk,i);
        
        PLOT=0;
        if PLOT
            
            figure;
            subplot(1,3,1);

            plot(xtmp,Etmp,'-b.',xtmp+GW(kk),Etmp,'-r.',...
                [0 max(x_data)+GW(end)],[Ei Ei],'k--',...
                x0,Ei,'mo',[x0 x0],[0 1],'k--',...
                x0,E_SLR(kk,i),'go',...
                [0 max(x_data)+GW(end)],[E_SLR(kk,i) E_SLR(kk,i)],'k--');
            
            axis([0 max(x_data)+GW(end) 0 1]);

            xlabel('Temp [C]')
            ylabel('Exceedance Prob.')

            subplot(1,3,2);
            
            plot(x_data,E_data,'-b.',x_data+GW(kk),E_data,'-r.',...
                [0.9*min(x_data) max(x_data)+GW(end)],[ri/TR ri/TR],'k--',...
                x0_data,E0,'mo',[x0_data x0_data],[0 1],'k--',...
                x0_data,E_SLR_DATA(kk,i),'go',...
                [0.9*min(x_data) max(x_data)+GW(end)],[E_SLR_DATA(kk,i) E_SLR_DATA(kk,i)],'k--',[x0 x0],[0 1],'g--');
            
            axis([0.9*min(x_data) max(x_data)+GW(end) 0 1]);
      
            subplot(1,3,3);
            
            plot(x1,E(x1,mu,sig,k),'-b.',x1,E(x1,mu+GW(kk),sig,k),'-r.',...
                [0.9*min(x_data) max(x_data)+GW(end)],[ri/TR ri/TR],'k--',...
                x_Tr(TR,mu,sig,k),E0,'mo',[x_Tr(TR,mu,sig,k) x_Tr(TR,mu,sig,k)],[0 1],'k--',...
                x_Tr(TR,mu,sig,k),E_SLR_GEV(kk,i),'go',...
                [0.9*min(x_data) max(x_data)+GW(end)],[E_SLR_GEV(kk,i) E_SLR_GEV(kk,i)],'k--');
            
            axis([0.9*min(x_data) max(x_data)+GW(end) 0 1]);
           
            xlabel('Temp [C]')
            ylabel('Exceedance Prob.')
            drawnow;
            stop

        end
        
    end
    
    ODDS_SLR(:,i)=E_SLR(:,i)./(1-E_SLR(:,i));
    ODDS_SLR_DATA(:,i)=E_SLR_DATA(:,i)./(1-E_SLR_DATA(:,i));
    ODDS_SLR_GEV(:,i)=E_SLR_GEV(:,i)./(1-E_SLR_GEV(:,i));
    
    ODDSi=Ei./(1-Ei);
    ODDS0=E0./(1-E0);
    
    color1=cmap(idx(i),:); % color based on clusters
    
    %color1='b';           % single color
    
    ALPHA=0.15;
    
    xp0=0.09;
    yp0=0.12;
    width=0.28;
    height=0.77;
    xsep=0.02;
    
    figure(44); % New line
    subplot('position',[xp0 yp0 width height]); hold on; box on;
    han1=plot(GW,ODDS_SLR(:,i)/ODDSi,'Color',color1); han1.Color(4)=ALPHA;  % set transparency
    set(han1,'LineWidth',2);
    title({'All data'},'Fontsize',14,'FontName',FontName);
    
    % single plot
    figure(444); hold on; box on;
    han1=plot(GW,ODDS_SLR(:,i)/ODDSi,'Color',color1); han1.Color(4)=ALPHA;  % set transparency
    set(han1,'LineWidth',2);

    figure(44);
    subplot('position',[xp0+(width+xsep) yp0 width height]); hold on; box on;
    han1=plot(GW,ODDS_SLR_DATA(:,i)/ODDS0,'Color',color1); han1.Color(4)=ALPHA;  % set transparency
    set(han1,'LineWidth',2);
    title({'Extreme events data'},'Fontsize',14,'FontName',FontName);
    
    subplot('position',[xp0+2*(width+xsep) yp0 width height]); hold on; box on;
    han1=plot(GW,ODDS_SLR_GEV(:,i)/ODDS0,'Color',color1); han1.Color(4)=ALPHA;  % set transparency
    set(han1,'LineWidth',2);
    title({'GEV model from extreme events data'},'Fontsize',14,'FontName',FontName);
    
    % find the amount of SLR to double odds
    
    xtmp=GW(idsel);
    
    ytmp=log2(ODDS_SLR(idsel,i)/ODDSi);      idn=~isnan(ytmp); p=polyfit(xtmp(idn),ytmp(idn),1);      % fit straight line  
    
    sigma_2x(i)=1/(nanmean(diff(ytmp)./diff(xtmp)));             % find sigma 2x
    sigma_2x_fit(i)=1/p(1);                                 % find sigma 2x using a regression
    
    ytmp=log2(ODDS_SLR_DATA(idsel,i)/ODDS0); idn=~isnan(ytmp); p=polyfit(xtmp(idn),ytmp(idn),1);      % fit straight line 
    
    sigma_2x_DATA(i)=1/(nanmean(diff(ytmp)./diff(xtmp)));       % find sigma 2x
    sigma_2x_fit_DATA(i)=1/p(1);                                % find sigma 2x using a regression
    
    ytmp=log2(ODDS_SLR_GEV(idsel,i)/ODDS0);  idn=~isnan(ytmp); p=polyfit(xtmp(idn),ytmp(idn),1);      % fit straight line  
    
    sigma_2x_GEV(i)=1/(nanmean(diff(ytmp)./diff(xtmp)));             % find sigma 2x
    sigma_2x_fit_GEV(i)=1/p(1);                                 % find sigma 2x using a regression

    if isnan(ODDSi)          % if the initial odds are undefined because the emperical probabiliy of the 50-year event cannot be interpolated... set sigma 2x to be a NaN
        sigma_2x(i)=NaN;
        sigma_2x_fit(i)=NaN;
    end

    PLOT2=1;
    if PLOT2
        
        han1=figure(4444); set(han1,'PaperpositionMode','auto','Position',[200 100 1200 700],'color','w'); hold on; box on;
        subplot('position',[0.09 0.12 0.41 0.80]); hold on; box on;
        han1=plot(GW,ODDS_SLR(:,i)/ODDSi,'Color',color1); han1.Color(4)=ALPHA;  % set transparency
        set(han1,'LineWidth',2);
        
        subplot('position',[0.58 0.12 0.41 0.80]); hold on; box on;
        han1=plot(GW/sigma_2x(i),log2(ODDS_SLR(:,i)/ODDSi),'Color',color1); han1.Color(4)=ALPHA;  % set transparency
        set(han1,'LineWidth',2);
        
        xx=GW/sigma_2x(i);
        yy=log2(ODDS_SLR(:,i)/ODDSi);
        
%         xx(isnan(yy))=[];
%         yy(isnan(yy))=[];
        
        xx(xx>8) = NaN; yy(xx>8) = NaN;
        yy(yy>8) = NaN; xx(yy>8) = NaN;

        for iii = 1:length(xx)
            if yy(iii)-xx(iii)>2 || yy(iii)-xx(iii)<-2
                nnnn = nnnn+1;
                break
            end
        end

%         if length(yy)>2
%             FAC_10(i)=interp1(xx,yy,10);
%         else
%             FAC_10(i)=NaN;
%         end
        
        figure(666); %hold on; box on;
        
        subplot(1,3,1);
        
        han1=plot(GW,ODDS_SLR(:,i)/ODDSi,'b',GW,2.^(GW/sigma_2x(i)),'r');
        set(han1,'LineWidth',2);
        set(gca,'Yscale','log','ytick',[1 2 4 8 16 32 64 128 256 512 1024],'YMinorTick','off','Xtick',0:0.1:0.5);
%         axis([0 SLR_max 1 1024]);
        axis([0 TR_max 1 1024]);
        
        subplot(1,3,2);
        
        han1=plot(GW,ODDS_SLR_DATA(:,i)/ODDS0,'b',GW,2.^(GW/sigma_2x_DATA(i)),'r');
        set(han1,'LineWidth',2);
        set(gca,'Yscale','log','ytick',[1 2 4 8 16 32 64 128 256 512 1024],'YMinorTick','off','Xtick',0:0.1:0.5);
%         axis([0 SLR_max 1 1024]);
        axis([0 TR_max 1 1024]);
        
        subplot(1,3,3);
        
        han1=plot(GW,ODDS_SLR_GEV(:,i)/ODDS0,'b',GW,2.^(GW/sigma_2x_GEV(i)),'r');
        set(han1,'LineWidth',2);
        set(gca,'Yscale','log','ytick',[1 2 4 8 16 32 64 128 256 512 1024],'YMinorTick','off','Xtick',0:0.1:0.5);
%         axis([0 SLR_max 1 1024]);
        axis([0 TR_max 1 1024]);
        
    end

    
end

% TR_max is the same as SLR_max
TR_max=max(GW); % TODO: ???? What this should be?
% TR_max = 10

% for i=1:3
%     subplot('position',[xp0+(i-1)*(width+xsep) yp0 width height]);
%     han1=plot(GW,2.^(GW/1),'--k',GW,2.^(GW/0.25),'--r',GW,2.^(GW/0.5),'--b',GW,2.^(GW/2),'--g');
%     set(han1,'LineWidth',3);
%     
% %     han1=plot([0 TR_max],[150 150],'--k');
%     
%     set(gca,'FontSize',14);
%     xlabel('Temperature rise [C]','FontSize',14);
%     if i==1;
%         ylabel('$\frac{E}{E_0}$','FontSize',30,'Interpreter','Latex'); 
%         set(gca,'Yscale','log','ytick',[1 2 4 8 16 32 64 128 256],'YMinorTick','off','Xtick',0:1:TR_max);
%     else
%         set(gca,'Yscale','log','ytick',[1 2 4 8 16 32 64 128 256],'yticklabel',[],'YMinorTick','off','Xtick',0:1:TR_max);
%     end
%     axis([0 TR_max 1 256]);
%     
%     text(1.00684527839933, 61.5049107158222,'0.25^o C','Rotation',85,'FontSize',14,'color','r');
%     text(2.54255952313542, 66.8957069204231,'0.5^o C','Rotation',75,'FontSize',14,'color','b');
%     text(4.80654761727367, 41.78037710968,'1^o C','Rotation',62,'FontSize',14,'color','k');
%     text(8.46696419375283, 23.4266560725073,'2^o C','Rotation',40,'FontSize',14,'color','g');
% 
% end


figure(44);
for i=1:3
    subplot('position',[xp0+(i-1)*(width+xsep) yp0 width height]);
%     han1=plot(GW,2.^(GW/0.05),'--k',GW,2.^(GW/0.01),'--r',GW,2.^(GW/0.25),'--g');
    plot(GW,2.^(GW/1),'--','LineWidth',3,'Color','k'); hold on
    plot(GW,2.^(GW/0.25),'--','LineWidth',3,'Color',[0.4940 0.1840 0.5560])
    plot(GW,2.^(GW/0.5),'--','LineWidth',3,'Color','b')
    plot(GW,2.^(GW/1.5),'--','LineWidth',3,'Color','r')

    
    %han1=plot([0 SLR_max],[150 150],'--k');
    
    set(gca,'FontSize',14,'FontName',FontName);
    xlabel('\mu_T [\circC]','FontSize',16,'FontName',FontName);
    if i==1
        ylabel('$\frac{\mathrm{Odds}}{\mathrm{Odds}_0}$','FontSize',26,'Interpreter','latex');
        set(gca,'Yscale','log','ytick',[1 2 4 8 16 32 64 128 256],'YMinorTick','off','Xtick',0:1:TR_max);
        % ????? Chenged it to the following
%         set(gca,'Yscale','log','ytick',[1 2 4 8],'YMinorTick','off','Xtick',0:0.1:0.5);
    else
        set(gca,'Yscale','log','ytick',[1 2 4 8 16 32 64 128 256],'yticklabel',[],'YMinorTick','off','Xtick',0:1:TR_max);
        % ????? Chenged it to the following
%         set(gca,'Yscale','log','ytick',[1 2 4 8],'yticklabel',[],'YMinorTick','off','Xtick',0:0.1:0.5);
    end
    axis([0 TR_max 1 256]);
%     axis([0 TR_max 1 8]);
    
%     text(0.0425595254415558,75.4266270583301,'1 C','Rotation',85,'FontSize',14,'color','r');
%     text(0.342261891989481,162.18128599867,'5 C','Rotation',57,'FontSize',14,'color','k');
%     text(0.401488106875193,3.60195568249456,'25 C','Rotation',18,'FontSize',14,'color','g');
    text(1.00684527839933, 61.5049107158222,'0.25 \circC','Rotation',80,'FontSize',14,'color',[0.4940 0.1840 0.5560],'FontName',FontName);
    text(2.48303569498516, 57.231344506145,'0.5 \circC','Rotation',75,'FontSize',14,'color','b','FontName',FontName);
    text(4.80654761727367, 41.78037710968,'1 \circC','Rotation',62,'FontSize',14,'color','k','FontName',FontName);
    text(6.59196436405182, 28.0478757570636,'1.5 \circC','Rotation',50,'FontSize',14,'color','r','FontName',FontName);

end

annotation(gcf,'textbox',[0.272500000000001 0.92 0.997500000000007 0.0633333333333341],'String',{'Normalized increase in the odds of exceedance (i.e., "doubling")'},'FitBoxToText','off','linestyle','none','FontWeight','bold','FontSize',18,'FontName',FontName);

annotation(gcf,'textbox',...
    [0.333833333333334 0.135000000000001 0.0286666666666665 0.0633333333333344],...
    'VerticalAlignment','middle',...
    'String',{'A'},...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',30,...
    'FontName','Nomada Incise',...
    'FitBoxToText','off',...
    'BackgroundColor',[0.8 0.8 0.8]);
annotation(gcf,'textbox',...
    [0.633833333333341 0.135000000000003 0.0286666666666665 0.0633333333333345],...
    'VerticalAlignment','middle',...
    'String','B',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',30,...
    'FontName','Nomada Incise',...
    'FitBoxToText','off',...
    'BackgroundColor',[0.8 0.8 0.8]);
annotation(gcf,'textbox',...
    [0.933833333333339 0.135000000000002 0.0286666666666664 0.0633333333333345],...
    'VerticalAlignment','middle',...
    'String','C',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',30,...
    'FontName','Nomada Incise',...
    'FitBoxToText','off',...
    'BackgroundColor',[0.8 0.8 0.8]);

% print(gcf,'figures/Figure4s_rate_of_increase_with_GTR_Appendix2.jpg','-djpeg','-r300');
exportgraphics(gcf,'figures/FigureA3.png','Resolution',300)


figure(444)

% han1=plot(GW,2.^(GW/0.05),'--k',GW,2.^(GW/0.01),'--r',GW,2.^(GW/0.25),'--g');
han1=plot(GW,2.^(GW/1),'--k',GW,2.^(GW/0.25),'--r',GW,2.^(GW/0.5),'--b',GW,2.^(GW/2),'--g');
set(han1,'LineWidth',3);

%han1=plot([0 SLR_max],[150 150],'--k');

set(gca,'FontSize',14);
xlabel('Global temperature rise [m]','FontSize',14);
ylabel('$\frac{\mathrm{Odds}}{\mathrm{Odds}_0}$','FontSize',30,'Interpreter','latex');
set(gca,'Yscale','log','ytick',[1 2 4 8 16 32 64 128 256],'YMinorTick','off','Xtick',0:1:TR_max);
axis([0 TR_max 1 256]);

% text(0.0425595254415558,75.4266270583301,'1 cm','Rotation',85,'FontSize',14,'color','r');
% text(0.342261891989481,162.18128599867,'5 cm','Rotation',57,'FontSize',14,'color','k');
% text(0.401488106875193,3.60195568249456,'25 cm','Rotation',18,'FontSize',14,'color','g');
plot(GW,2.^(GW/1),'--','LineWidth',3,'Color','k'); hold on
plot(GW,2.^(GW/0.25),'--','LineWidth',3,'Color',[0.4940 0.1840 0.5560])
plot(GW,2.^(GW/0.5),'--','LineWidth',3,'Color','b')
plot(GW,2.^(GW/1.5),'--','LineWidth',3,'Color','r')

annotation(gcf,'textbox',[0.0384615384615385 0.941666666666667 0.999999999999999 0.0633333333333341],'String',{'Normalized increase in odds of exceedance (i.e., "doubling")'},'FitBoxToText','off','linestyle','none','FontWeight','bold','FontSize',16);

% exportgraphics(gcf,'figures/FigureA3.png','Resolution',300)

figure(4444);
subplot('position',[0.09 0.12 0.41 0.80]); hold on; box on;

% han1=plot(GW,2.^(GW/0.05),'--k',GW,2.^(GW/0.01),'--r',GW,2.^(GW/0.25),'--g');
plot(GW,2.^(GW/1),'--','LineWidth',3,'Color','k'); hold on
plot(GW,2.^(GW/0.25),'--','LineWidth',3,'Color',[0.4940 0.1840 0.5560])
plot(GW,2.^(GW/0.5),'--','LineWidth',3,'Color','b')
plot(GW,2.^(GW/1.5),'--','LineWidth',3,'Color','r')

set(gca,'FontSize',18,'FontName',FontName);
xlabel('\mu_T [\circC]','FontSize',18,'FontName',FontName);
ylabel('$\frac{\mathrm{Odds}}{\mathrm{Odds}_0}$','FontSize',28,'Interpreter','latex');
set(gca,'Yscale','log','ytick',[1 2 4 8 16 32 64 128 256],'YMinorTick','off','Xtick',0:1:TR_max);
axis([0 TR_max 1 256]);

text(1.2214212446678, 74.6834240517126,'0.25 \circC','Rotation',80,'FontSize',17,'color',[0.4940 0.1840 0.5560],'FontName',FontName);
text(2.70331748132783, 69.7233636218025,'0.5 \circC','Rotation',72,'FontSize',17,'color','b','FontName',FontName);
text(5.23453807927729, 51.064070450407,'1 \circC','Rotation',57,'FontSize',17,'color','k','FontName',FontName);
text(6.80900867392377, 28.8828549473691,'1.5 \circC','Rotation',45,'FontSize',17,'color','r','FontName',FontName);

subplot('position',[0.58 0.12 0.41 0.80]); hold on; box on;
han1=plot([0 45],[0 45],'--k'); set(han1,'LineWidth',3);
han1=plot([0 45],[0+2 45+2],'--k'); set(han1,'LineWidth',1);
han1=plot([2 45],[0 45-2],'--k'); set(han1,'LineWidth',1);
% axis([0 15 0 15]);
axis([0 8 0 8]);

annotation(gcf,'textbox',...
    [0.236000000000004 0.918809523809525 0.607333333333336 0.0633333333333342],...
    'String','Normalized increase in the odds of exceedance (i.e., "doubling")',...
    'LineStyle','none',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName',FontName,...
    'FitBoxToText','off');
annotation(gcf,'textbox',...
    [0.464166666666667 0.127571428571429 0.0316666666666678 0.0552857142857144],...
    'VerticalAlignment','middle',...
    'String','A',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',30,...
    'FontName','Nomada Incise',...
    'FitBoxToText','off',...
    'BackgroundColor',[0.8 0.8 0.8]);
annotation(gcf,'textbox',...
    [0.954166666666673 0.127571428571428 0.031666666666668 0.0552857142857143],...
    'VerticalAlignment','middle',...
    'String','B',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',30,...
    'FontName','Nomada Incise',...
    'FitBoxToText','off',...
    'BackgroundColor',[0.8 0.8 0.8]);
set(gca,'ytick',0:2:16,'xtick',0:2:16,'FontSize',18,'FontName',FontName);

xlabel('$\frac{\mu_\mathrm{T}}{\tilde{\sigma}}$','FontSize',28,'interpreter','latex');
ylabel('$\log2\left(\frac{\mathrm{Odds}}{\mathrm{Odds}_0}\right)$','FontSize',22,'Interpreter','latex');

annotation(gcf,'textarrow',[0.794687500000002 0.695],...
    [0.21984126984127 0.331428571428571],'String',{'$\frac{\mathrm{Odds}}{\mathrm{Odds}_0}=2^{\frac{\mu_\mathrm{T}}{\tilde{\sigma}}}$'},'FontSize',30,'interpreter','latex');

% print(gcf,'figures/Figure4_rate_of_increase_with_GTR_ODDS_ALTERNATE.jpg','-djpeg','-r300');
exportgraphics(gcf,'figures/Figure5.png','Resolution',300)

% FAC_10(isnan(FAC_10))=[];
% sum(FAC_10>8 & FAC_10<12)./length(FAC_10)*100


%% plot sigma 2x (Figure 5)

han1=figure(31); set(han1,'PaperpositionMode','auto','Position',[200 100 1000 400],'color','w'); hold on; box on;

id=~isnan(sigma_2x);

subplot('position',[0.056 0.11 0.542 0.8]); hold on; box on;
geoshow('landareas.shp', 'FaceColor',[0.6 0.6 0.6],'EdgeColor',[0.6 0.6 0.6]); 
axis equal; axis tight;

scatter(lon(id),lat(id),30*ones(size(lat(id))),sigma_2x(id),'filled'); 
colormap(gca,flipud(jet)); 
caxis([0.2 1.8]); 
colorbar;

axis([-180 -60 0 80]);
set(gca,'FontSize',14,'FontName',FontName);
title('Spatial variation in  ','FontName',FontName,'FontSize',18)

% id_west_coast=find(lon<-110)
% id_east_coast=find(lon>-110)
% 
% nanmean(sigma_2x(id_west_coast))
% nanmean(sigma_2x(id_east_coast))

subplot('position',[0.69 0.11 0.3 0.8]); hold on; box on;
% xe=(0:0.015:0.25)';
xe=(0:0.15:4)';

ybar=NaN(length(xe),4);

% cmap2=cmap([3 4 2 1],:);
cmap2=cmap([2 3 1 4],:);


for i=1:length(xe)-1
   ybar(i,1)=sum(idx==2 & sigma_2x>=xe(i) & sigma_2x<xe(i+1));
   ybar(i,2)=sum(idx==3 & sigma_2x>=xe(i) & sigma_2x<xe(i+1));
   ybar(i,3)=sum(idx==1 & sigma_2x>=xe(i) & sigma_2x<xe(i+1));
   ybar(i,4)=sum(idx==4 & sigma_2x>=xe(i) & sigma_2x<xe(i+1));
end

ybar(end,1)=sum(idx==2 & sigma_2x>xe(end));
ybar(end,2)=sum(idx==3 & sigma_2x>xe(end));
ybar(end,3)=sum(idx==1 & sigma_2x>xe(end));
ybar(end,4)=sum(idx==4 & sigma_2x>xe(end));

bh=bar(xe,ybar,'stacked'); 
for i=1:4
    bh(i).FaceColor = cmap2(i,:);
end

text(0.77, 1100,['$\tilde{\sigma}$ (mean)    = ~' num2str(nanmean(sigma_2x),2),' $^\circ$C'],'FontSize',18,'interpreter','latex');
text(0.77, 1000,['$\tilde{\sigma}$ (median) = ~' num2str(nanmedian(sigma_2x),2),' $^\circ$C'],'FontSize',18,'interpreter','latex');
text(0.77, 900,['$\tilde{\sigma}$ (s.d.)       = ~' num2str(nanstd(sigma_2x),2),' $^\circ$C'],'FontSize',18,'interpreter','latex');

ylabel('Count (No. of stations)','FontSize',14,'FontName',FontName);
title('Histogram plot of  ','FontName',FontName,'FontSize',18);

% set(gca,'FontSize',16,'XLim',[0 0.26],'XTick',[0:0.05:0.25]);
set(gca,'FontSize',14,'XLim',[0 2],'XTick',[0:0.5:2],'FontName',FontName);
xlabel(append('$\tilde{\sigma}$', ' [$^\circ$C]'),'FontSize',18,'FontName',FontName,'Interpreter','latex');

annotation(gcf,'textbox',[0.55 0.9225 0.088 0.075],'String',{'  [\circC]'},'FitBoxToText','off','LineStyle','none','FontSize',18,'FontName',FontName);
annotation(gcf,'textbox',[0.542 0.9175 0.088 0.075],...
    'String',{'$\tilde{\sigma}$'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',18,...
    'FitBoxToText','off');
annotation(gcf,'textbox',[0.387 0.905 0.0909999999999998 0.075],...
    'String',{'$\tilde{\sigma}$'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',20,...
    'FitBoxToText','off');
annotation(gcf,'textbox',...
    [0.907000000000001 0.9 0.0910000000000007 0.075],...
    'String',{'$\tilde{\sigma}$'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',20,...
    'FitBoxToText','off');


annotation(gcf,'textbox',...
    [0.0651666666666672 0.7975 0.0378333333333328 0.0978571428571433],...
    'VerticalAlignment','middle',...
    'String','A',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',30,...
    'FontName','Nomada Incise',...
    'FitBoxToText','off',...
    'BackgroundColor',[0.8 0.8 0.8]);
annotation(gcf,'textbox',...
    [0.696166666666668 0.7975 0.0378333333333328 0.0978571428571433],...
    'VerticalAlignment','middle',...
    'String','B',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',30,...
    'FontName','Nomada Incise',...
    'FitBoxToText','off',...
    'BackgroundColor',[0.8 0.8 0.8]);

% print(gcf,'figures/Figure5_GTR_SCALE.jpg','-djpeg','-r300');
exportgraphics(gcf,'figures/Figure6.png','Resolution',300)


%% plot finc

figure(45)

sig_tilde1=repmat(GW,1,len)./log2(E_SLR./repmat(EI,length(GW),1));
sig_tilde2=repmat(GW,1,len)./log2(E_SLR_DATA./E0);
sig_tilde3=repmat(GW,1,len)./log2(E_SLR_GEV./E0);

finc1=nanmean(log2(E_SLR./repmat(EI,length(GW),1)),2);
finc2=nanmean(log2(E_SLR_DATA./E0),2);
finc3=nanmean(log2(E_SLR_GEV./E0),2);

plot(GW,finc1,'r',GW,finc2,'g',GW,finc3,'b');

id=find(GW>0.1,1,'first');

sig_tilde_bar1=nanmean(sig_tilde1(1:id,:));
sig_tilde_bar2=nanmean(sig_tilde2(1:id,:));
sig_tilde_bar3=nanmean(sig_tilde3(1:id,:));

han1=figure(444); set(han1,'PaperpositionMode','auto','Position',[200 100 1200 500],'color','w');

x0=0.08;
xsep=0.025;
y0=0.08;
width=0.28;
height=0.9;

AXMAX=0.6

subplot('Position',[x0 y0 width height]); hold on; box on;
scatter([S(ID).sigma],sig_tilde_bar1,30*ones(len,1),[S(ID).k],'filled'); colormap(jet(64)); 
plot([0 AXMAX],[0 AXMAX],'k--'); axis equal;
axis([0 AXMAX 0 AXMAX]);
set(gca,'FontSize',14);
xlabel('$\sigma$','FontSize',30,'Interpreter','Latex');
ylabel('$\tilde{\sigma}$','FontSize',30,'Interpreter','Latex');
 title('for all data');

subplot('Position',[x0+(xsep+width) y0 width height]); hold on; box on;
scatter([S(ID).sigma],sig_tilde_bar2,30*ones(len,1),[S(ID).k],'filled'); colormap(jet(64)); 
plot([0 AXMAX],[0 AXMAX],'k--'); axis equal;
axis([0 AXMAX 0 AXMAX]);
set(gca,'YTickLabel',[],'FontSize',14);
xlabel('$\sigma$','FontSize',30,'Interpreter','Latex');
title('for extreme events data');

subplot('Position',[x0+2*(xsep+width) y0 width height]); hold on; box on;
scatter([S(ID).sigma],sig_tilde_bar3,30*ones(len,1),[S(ID).k],'filled'); colormap(jet(64)); colorbar
plot([0 AXMAX],[0 AXMAX],'k--'); axis equal;
axis([0 AXMAX 0 AXMAX]);
set(gca,'YTickLabel',[],'FontSize',14);
xlabel('$\sigma$','FontSize',30,'Interpreter','Latex');
title('for GEV model of extremes');

ss1=[S(ID).sigma];
kk1=[S(ID).k];

id=abs(kk1)<0.01;
plot(ss1(id),sig_tilde_bar3(id),'ko');

annotation(gcf,'textbox',[0.09 0.78 0.0236666666666665 0.07],'String',{'A'},'FitBoxToText','off','FontSize',16,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle');
annotation(gcf,'textbox',[0.09+1*(width+xsep) 0.78 0.0236666666666665 0.07],'String',{'B'},'FitBoxToText','off','FontSize',16,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle');
annotation(gcf,'textbox',[0.09+2*(width+xsep) 0.78 0.0236666666666665 0.07],'String',{'C'},'FitBoxToText','off','FontSize',16,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle');


%% Exceedance Plot with Time

% quick check of SLR curves

% t_Temp=datenum(2000,1,1):365.25:datenum(2100,1,1);
% % TEMP
% 
% 
% han1=figure(5); set(han1,'PaperpositionMode','auto','Position',[200 100 1200 600],'color','w');
% 
% TIME=(datenum(2000,1,1):100:datenum(2300,1,1))';
%     
% E_TIME=NaN(length(TIME),len);
% 
% E_TIME_DATA=NaN(length(TIME),len);
% 
% E_TIME_GEV=NaN(length(TIME),len);
%     
% TR=50;
% 
% x1=(0:0.001:2*max(x_data))';
% 
% % for i=1:1:len
% for i=1:10:len This will not go over all the stations
%     
%     id=ID(i); % get the specific station
%     
%     E0=S(id).ri/TR; % get the exceedance probability of the TR-year event
%     
%     % GEV
%     mu=S(id).mu;
%     sig=S(id).sigma;
%     k=S(id).k;
% 
% %     GW=interp1(S(id).t_SLR,S(id).SLR_med(:,SCENARIO),TIME);
%     GW=interp1(ti,Tmean,TIME);
% 
%     %plot(S(id).t_SLR,S(id).SLR_med(:,SCENARIO),'ro',TIME,SLR,'-b.'); datetick('x');
%     
%     % empirical exceedace probability of all data
%     [Etmp,ia,~] = unique(S(id).Exceedance,'last');
%     xtmp=S(id).x(ia);
% 
%     x0=x_Tr(TR,mu,sig,k); % find the water level corresponding to the 50-year event
%     
%     Ei=interp1(xtmp,Etmp,x0,'linear');
% 
%     % empirical exceedace probability of the data of extreme events
%     len_data=length(S(id).h_max);
%     E_data=linspace(1/(len_data+1),len_data/(len_data+1),len_data)';
%     x_data=S(id).h_max;
% 
%     x0_data=interp1(E_data,x_data,E0,'linear','extrap');
%     %Ei_data=interp1(x_data,E_data,x0,'linear');
%     
%     [tmp_x,~,~] = unique(x_data,'last');
%     tmp_y=tmp_x;
%     for jj=1:length(tmp_x)
%         ind_tmp=find(x_data==tmp_x(jj));
%         tmp_y(jj)=mean(E_data(ind_tmp));
%     end
%     
%     for kk=1:length(GW)       
%         % find empirical exceedance probability increase of the all values
%         E_TIME(kk,i)=interp1(xtmp+GW(kk),Etmp,x0);
% 
%         % find empirical exceedance probability increase of the extremes
% %         E_TIME_DATA(kk,i)=interp1(x_data+GW(kk),E_data,x0_data);
%                
%         E_TIME_DATA(kk,i)=interp1(tmp_x+GW(kk),tmp_y,x0_data);
% 
%         % find GEV probability increase
%         E_TIME_GEV(kk,i)=E(x0,mu+GW(kk),sig,k);
%  
%         
%         PLOT=0;
%         
%         if E_TIME_DATA(kk,i)<0
%             PLOT=1;
%         end
%         
%         if PLOT
%             
%             figure;
%             subplot(1,3,1);
% 
%             plot(xtmp,Etmp,'-b.',xtmp+GW(kk),Etmp,'-r.',...
%                 [0 max(x_data)+GW(end)],[Ei Ei],'k--',...
%                 x0,Ei,'mo',[x0 x0],[0 1],'k--',...
%                 x0,E_TIME(kk,i),'go',...
%                 [0 max(x_data)+GW(end)],[E_TIME(kk,i) E_TIME(kk,i)],'k--');
%             
%             axis([0 max(x_data)+GW(end) 0 1]);
% 
%             subplot(1,3,2);
%             
%             plot(x_data,E_data,'-b.',x_data+GW(kk),E_data,'-r.',...
%                 [0.9*min(x_data) max(x_data)+GW(end)],[S(id).ri/TR S(id).ri/TR],'k--',...
%                 x0_data,E0,'mo',[x0_data x0_data],[0 1],'k--',...
%                 x0_data,E_TIME_DATA(kk,i),'go',...
%                 [0.9*min(x_data) max(x_data)+GW(end)],[E_TIME_DATA(kk,i) E_TIME_DATA(kk,i)],'k--',[x0 x0],[0 1],'g--');
%             
%             axis([0.9*min(x_data) max(x_data)+GW(end) 0 1]);
%       
%             subplot(1,3,3);
%             
%             plot(x1,E(x1,mu,sig,k),'-b.',x1,E(x1,mu+GW(kk),sig,k),'-r.',...
%                 [0.9*min(x_data) max(x_data)+GW(end)],[S(id).ri/TR S(id).ri/TR],'k--',...
%                 x_Tr(TR,mu,sig,k),E0,'mo',[x_Tr(TR,mu,sig,k) x_Tr(TR,mu,sig,k)],[0 1],'k--',...
%                 x_Tr(TR,mu,sig,k),E_TIME_GEV(kk,i),'go',...
%                 [0.9*min(x_data) max(x_data)+GW(end)],[E_TIME_GEV(kk,i) E_TIME_GEV(kk,i)],'k--');
%             
%             axis([0.9*min(x_data) max(x_data)+GW(end) 0 1]);
%            
%             drawnow;
%             stop
% 
%         end
%         
%     end
% 
%     color1=cmap(idx(i),:); % color based on clusters
% 
% %     ????
% %     if ismember(i,id_neg)
% %         color1='k';           % single color
% %     end
%     
%     ALPHA=0.15;
%     
%     xp0=0.09;
%     yp0=0.12;
%     width=0.28;
%     height=0.77;
%     xsep=0.02;
%     
%     subplot('position',[xp0 yp0 width height]); hold on; box on;
%     han1=plot(TIME,E_TIME(:,i)/Ei,'Color',color1); han1.Color(4)=ALPHA;  % set transparency
%     set(han1,'LineWidth',2);
%     title({'for all data'});
%     datetick('x');
%     
%     subplot('position',[xp0+(width+xsep) yp0 width height]); hold on; box on;
%     han1=plot(TIME,E_TIME_DATA(:,i)/E0,'Color',color1); han1.Color(4)=ALPHA;  % set transparency
%     set(han1,'LineWidth',2);
%     title({'for extreme events data'});
%     datetick('x');
%     
%     subplot('position',[xp0+2*(width+xsep) yp0 width height]); hold on; box on;
%     han1=plot(TIME,E_TIME_GEV(:,i)/E0,'Color',color1); han1.Color(4)=ALPHA;  % set transparency
%     set(han1,'LineWidth',2);
%     title({'for GEV model of extremes'})
%     datetick('x');
% 
% end
% 
% TIME_min=datenum(2000,1,1);
% TIME_max=datenum(2300,1,1);
% 
% annotation(gcf,'textbox',[0.268500000000001 0.941666666666667 0.607333333333333 0.0633333333333341],'String',{'Normalized exceedance probability increase (i.e., "doubling")'},'FitBoxToText','off','linestyle','none','FontWeight','bold','FontSize',16);
% 
% annotation(gcf,'textbox',[0.340500000000001 0.825000000000001 0.0236666666666666 0.0583333333333342],'String',{'A'},'FontWeight','bold','FontSize',16,'FitBoxToText','off','BackgroundColor',[1 1 1],'HorizontalAlignment','center','VerticalAlignment','middle');
% annotation(gcf,'textbox',[0.340500000000001+1*(width+xsep) 0.825000000000001 0.0236666666666666 0.0583333333333342],'String',{'B'},'FontWeight','bold','FontSize',16,'FitBoxToText','off','BackgroundColor',[1 1 1],'HorizontalAlignment','center','VerticalAlignment','middle');
% annotation(gcf,'textbox',[0.340500000000001+2*(width+xsep) 0.825000000000001 0.0236666666666666 0.0583333333333342],'String',{'C'},'FontWeight','bold','FontSize',16,'FitBoxToText','off','BackgroundColor',[1 1 1],'HorizontalAlignment','center','VerticalAlignment','middle');
% 
% for i=1:3    
% 
%     subplot('position',[xp0+(i-1)*(width+xsep) yp0 width height]);
%     
%     text(747180.471633033, 93.0160170319062,'5 yr','Rotation',85,'FontSize',14,'color','r');
%     text(758920.317284235, 71.5733621237796,'10 yrs','Rotation',80,'FontSize',14,'color','k');
%     text(791866.558016729, 70.0335658386506,'25 yrs','Rotation',65,'FontSize',14,'color','b');
% %     text(770210.538767274,14.6392415674771,'25 yrs','Rotation',55,'FontSize',14,'color','b','BackgroundColor',[1 1 1]);
%     
% %     han1=plot(TIME,2.^((TIME-TIME_min)/(5*365.25)),'--k',TIME,2.^((TIME-TIME_min)/365.25),'--r',TIME,2.^((TIME-TIME_min)/(25*365.25)),'--g');
%     han1=plot(TIME,2.^((TIME-TIME_min)/(10*365.25)),'--k',TIME,2.^((TIME-TIME_min)/(5*365.25)),'--r',TIME,2.^((TIME-TIME_min)/(25*365.25)),'--b');
%     set(han1,'LineWidth',3);
%     
%     han1=plot([TIME_min TIME_max],[150 150],'--k');
%     
%     set(gca,'FontSize',14);
%     if i==1;
%         ylabel('$\frac{E}{E_0}$','FontSize',30,'Interpreter','Latex'); 
%         set(gca,'Yscale','log','ytick',[1 2 4 8 16 32 64 128 256],'YMinorTick','off');
%     else
%         set(gca,'Yscale','log','ytick',[1 2 4 8 16 32 64 128 256],'yticklabel',[],'YMinorTick','off');
%     end
%     
%     if i>1;
%         set(gca,'xtick',datenum(2000:50:2300,1,1),...
%             'XTickLabel',{'' '2050' '2100' '2150' '2200' '2250' '2300'})
%     end
%     axis([TIME_min TIME_max 1 256]);
%     datetick
% end

% print(gcf,'figures/Figure5_rate_of_increase_with_time.jpg','-djpeg','-r300');

%% Exceedance Plot with Time (Appendix 3)

% quick check of SLR curves

han1=figure(6); set(han1,'PaperpositionMode','auto','Position',[200 100 1600 1000],'color','w');

TIME_min=datenum(2000,1,1);
% TIME_max=datenum(2300,1,1);
% The new changes
TIME_max=datenum(2100,1,1);
TEMP_PRJ=[Tmin, Tmean, Tmax]; % ????

for SCENARIO=1:3
    
    for i=1:1:len
        
%         id=ID(i); % get the specific station
%         [SLR_min(i),tid]=min(S(id).SLR_med(:,SCENARIO));
        [TEMP_min(i),tid]=min(TEMP_PRJ(:,SCENARIO));
        
        %plot(S(id).t_SLR,S(id).SLR_med(:,SCENARIO),'b-',S(id).t_SLR(tid),S(id).SLR_med(tid,SCENARIO),'ro'); datetick('x');
    end
    
    id_neg=find(TEMP_min<-0.01);
    
    % S(ID(id_neg)).name
    
%     TIME=(datenum(2000,1,1):100:datenum(2300,1,1))';
%   The new changes:
%     TIME=(datenum(2000,1,1):100:datenum(2150,1,1))'; % Global
    TIME=(T_proj_local.time_start:100:datenum(2100,1,1))'; % Local
    
    E_TIME=NaN(length(TIME),len);
    
    E_TIME_DATA=NaN(length(TIME),len);
    
    E_TIME_GEV=NaN(length(TIME),len);
    
    TR=50;
    
    x1=(0:0.001:2*max(x_data))';
    
%     for i=1:1:len
    for i=1:1:len %  We are not going over all the stations 
        [SCENARIO i]
        
        id=ID(i); % get the specific station
        
        E0=S(id).ri/TR; % get the exceedance probability of the TR-year event
        
        % GEV
        mu=S(id).mu;
        sig=S(id).sigma;
        k=S(id).k;
        
%         GW=interp1(S(id).t_SLR,S(id).SLR_med(:,SCENARIO),TIME);
%       NWE change:
        GW=interp1(S(id).ti,S(id).TR_scenario(:,SCENARIO), TIME); %local 
%         GW=interp1(ti,TEMP_PRJ(:,SCENARIO),TIME); % Global
        
        %plot(S(id).t_SLR,S(id).SLR_med(:,SCENARIO),'ro',TIME,SLR,'-b.'); datetick('x');
        
        % empirical exceedace probability of all data
        [Etmp,ia,~] = unique(S(id).Exceedance,'last');
        xtmp=S(id).x(ia);
        
        x0=x_Tr(TR,mu,sig,k); % find the water level corresponding to the 50-year event
        
        Ei=interp1(xtmp,Etmp,x0,'linear');
        
        % empirical exceedace probability of the data of extreme events
        len_data=length(S(id).h_max);
        E_data=linspace(1/(len_data+1),len_data/(len_data+1),len_data)';
        x_data=S(id).h_max;
        
        x0_data=interp1(E_data,x_data,E0,'linear','extrap');
        %x0_data=interp1(E_data,x_data,E0);
        %Ei_data=interp1(x_data,E_data,x0,'linear');
        [tmp_x,~,~] = unique(x_data,'last');
        tmp_y=tmp_x;
        for jj=1:length(tmp_x)
            ind_tmp=find(x_data==tmp_x(jj));
            tmp_y(jj)=mean(E_data(ind_tmp));
        end
        
        for kk=1:length(GW)
            
            % find empirical exceedance probability increase of the all values
            E_TIME(kk,i)=interp1(xtmp+GW(kk),Etmp,x0);
            
            % find empirical exceedance probability increase of the extremes
%             E_TIME_DATA(kk,i)=interp1(x_data+GW(kk),E_data,x0_data);
            E_TIME_DATA(kk,i)=interp1(tmp_x+GW(kk),tmp_y,x0_data);
            
            % find GEV probability increase
            E_TIME_GEV(kk,i)=E(x0,mu+GW(kk),sig,k);
            
            PLOT=0;
            
            if E_TIME_DATA(kk,i)<0
                PLOT=1;
            end
            
            if PLOT
                
                figure;
                subplot(1,3,1);
                
                plot(xtmp,Etmp,'-b.',xtmp+GW(kk),Etmp,'-r.',...
                    [0 max(x_data)+GW(end)],[Ei Ei],'k--',...
                    x0,Ei,'mo',[x0 x0],[0 1],'k--',...
                    x0,E_TIME(kk,i),'go',...
                    [0 max(x_data)+GW(end)],[E_TIME(kk,i) E_TIME(kk,i)],'k--');
                
                axis([0 max(x_data)+GW(end) 0 1]);
                
                subplot(1,3,2);
                
                plot(x_data,E_data,'-b.',x_data+GW(kk),E_data,'-r.',...
                    [0.9*min(x_data) max(x_data)+GW(end)],[S(id).ri/TR S(id).ri/TR],'k--',...
                    x0_data,E0,'mo',[x0_data x0_data],[0 1],'k--',...
                    x0_data,E_TIME_DATA(kk,i),'go',...
                    [0.9*min(x_data) max(x_data)+GW(end)],[E_TIME_DATA(kk,i) E_TIME_DATA(kk,i)],'k--',[x0 x0],[0 1],'g--');
                
                axis([0.9*min(x_data) max(x_data)+GW(end) 0 1]);
                
                subplot(1,3,3);
                
                plot(x1,E(x1,mu,sig,k),'-b.',x1,E(x1,mu+GW(kk),sig,k),'-r.',...
                    [0.9*min(x_data) max(x_data)+GW(end)],[S(id).ri/TR S(id).ri/TR],'k--',...
                    x_Tr(TR,mu,sig,k),E0,'mo',[x_Tr(TR,mu,sig,k) x_Tr(TR,mu,sig,k)],[0 1],'k--',...
                    x_Tr(TR,mu,sig,k),E_TIME_GEV(kk,i),'go',...
                    [0.9*min(x_data) max(x_data)+GW(end)],[E_TIME_GEV(kk,i) E_TIME_GEV(kk,i)],'k--');
                
                axis([0.9*min(x_data) max(x_data)+GW(end) 0 1]);
                
                drawnow;
                stop
                
            end
            
        end
        
        color1=cmap(idx(i),:); % color based on clusters
        
        if ismember(i,id_neg)
            color1='k';           % single color
        end
        
        ALPHA=0.15;
        
        xp0=0.045;
        yp0=0.66;
        width=0.205;
        height=0.27;
        xsep1=0.05;
        xsep2=0.025;
        ysep=0.03;
        
        
        
        subplot('position',[xp0 yp0-(SCENARIO-1)*(height+ysep) width height]); hold on; box on;
        han1=plot(TIME,GW,'Color',color1); han1.Color(4)=ALPHA;  % set transparency
        set(han1,'LineWidth',2);

        
        subplot('position',[xp0+xsep1+1*(width+xsep2) yp0-(SCENARIO-1)*(height+ysep) width height]); hold on; box on;
        han1=plot(TIME,E_TIME(:,i)/Ei,'Color',color1); han1.Color(4)=ALPHA;  % set transparency
        set(han1,'LineWidth',2);
        if SCENARIO==1
            title({'All data'},'FontSize',18,'FontName',FontName);
        end
        set(gca,'FontSize',14,'FontName',FontName);
        datetick('x');

        subplot('position',[xp0+xsep1+2*(width+xsep2) yp0-(SCENARIO-1)*(height+ysep) width height]); hold on; box on;
        han1=plot(TIME,E_TIME_DATA(:,i)/E0,'Color',color1); han1.Color(4)=ALPHA;  % set transparency
        set(han1,'LineWidth',2);
        if SCENARIO==1
            title({'Extreme events data'},'FontSize',18,'FontName',FontName);
        end
        set(gca,'FontSize',14,'FontName',FontName);
        datetick('x');
        
        subplot('position',[xp0+xsep1+3*(width+xsep2) yp0-(SCENARIO-1)*(height+ysep) width height]); hold on; box on;
        han1=plot(TIME,E_TIME_GEV(:,i)/E0,'Color',color1); han1.Color(4)=ALPHA;  % set transparency
        set(han1,'LineWidth',2);
        if SCENARIO==1
            title({'GEV model from extreme events data'},'FontSize',18,'FontName',FontName)
        end
        set(gca,'FontSize',14,'FontName',FontName);
        datetick('x');
        
    end
    
    % Average GW
    GW=interp1(ti,TEMP_PRJ(:,SCENARIO),TIME);
    subplot('position',[xp0 yp0-(SCENARIO-1)*(height+ysep) width height]); hold on; box on;
%     han1=plot(TIME,GW,'Color','b'); %han1.Color(4)=ALPHA;  % set transparency
%     set(han1,'LineWidth',2);
    
    
%     if SCENARIO==1
%         title({'ensemble mean GTR'});
%         ylabel('Best temp. scenario [C]','FontSize',14);
%     elseif SCENARIO==2
%         ylabel('Average temp. scenario [C]','FontSize',14);
%     elseif SCENARIO==3
%         ylabel('Worst temp. scenario [C]','FontSize',14);
%     end

    if SCENARIO==1
        title({'Local temperature projections'},'FontSize',18,'FontName',FontName)
        ylabel('Low-end RCP 8.5 TR [\circC]','FontSize',14);
    elseif SCENARIO==2
        ylabel('Intermediate RCP 8.5 TR [\circC]','FontSize',14);
    elseif SCENARIO==3
        ylabel('High-end RCP 8.5 TR [\circC]','FontSize',14);
    end

%     axis([TIME_min TIME_max -0.5 15]);
    axis([TIME_min TIME_max -0.5 10.5]);
    
    
    set(gca,'FontSize',14,'FontName',FontName);
    datetick('x');
    xlabel('Year','FontName',FontName);
    if SCENARIO<3
        set(gca,'xticklabel',[]);
        xlabel('');
    end

    if i == i(end)
        if SCENARIO == 1
            AVG_Global = Tmin(1:1212,:);
            AVG_Local = mean(T_low_stations(:,1:101),1);
        elseif SCENARIO == 2
            AVG_Global = Tmean(1:1212,:);
            AVG_Local = mean(T_mean_stations(:,1:101),1);
        elseif SCENARIO == 3
            AVG_Global = Tmax(1:1212,:);
            AVG_Local = mean(T_high_stations(:,1:101),1);
        end

        AVG_Global_int = interp1(1:length(AVG_Global),AVG_Global',1:length(AVG_Global)/length(TIME):length(AVG_Global));
        AVG_Local_int = interp1(1:length(AVG_Local),AVG_Local,linspace(1,length(AVG_Local),length(TIME)));

        han11 = plot(TIME,AVG_Global_int,'--k','LineWidth',2);
        han22 = plot(TIME,AVG_Local_int,'--b','LineWidth',2);

        if SCENARIO == 1
            legend([han11 han22],'Global average','North America average','FontSize',14,'FontName',FontName);
        end
    end


    TIME_min=datenum(2000,1,1);
%     TIME_max=datenum(2300,1,1);
    TIME_max=datenum(2100,1,1);
    
    for i=1:3
        
        subplot('position',[xp0+xsep1+i*(width+xsep2) yp0-(SCENARIO-1)*(height+ysep) width height]); hold on; box on;
        
%         text(737986.09325084, 73.303021521361,'5 yr','Rotation',85,'FontSize',14,'color','r');
%         text(746077.578529558, 48.9760389940741 ,'10 yrs','Rotation',70,'FontSize',14,'color','k');
%         text(773144.572316924, 16.2514969684399,'25 yrs','Rotation',50,'FontSize',14,'color','b');%,'BackgroundColor',[1 1 1]

%         New changes
        text(731306.861571456, 82.6691307776895,'1 year','Rotation',88,'FontSize',14,'color','m','FontName',FontName);
        text(739435.75474142, 59.3537289020805 ,'5 years','Rotation',67,'FontSize',14,'color','b','FontName',FontName);
        text(746790.844100608, 33.7572473193161 ,'10 years','Rotation',50,'FontSize',14,'color','k','FontName',FontName);
        text(755415.747120308, 8.80139136611196 ,'25 years','Rotation',25,'FontSize',14,'color',[0.4940 0.1840 0.5560],'FontName',FontName);
        text(757698.809620308, 3.64269028542452 ,'50 years','Rotation',15,'FontSize',14,'color','r','FontName',FontName);

%         han1=plot(TIME,2.^((TIME-TIME_min)/(5*365.25)),'--k',TIME,2.^((TIME-TIME_min)/365.25),'--r',TIME,2.^((TIME-TIME_min)/(25*365.25)),'--g');

        plot(TIME,2.^((TIME-TIME_min)/(5*365.25)),'--','LineWidth',3,'Color','b'); hold on
        plot(TIME,2.^((TIME-TIME_min)/365.25),'--','LineWidth',3,'Color','m')
        plot(TIME,2.^((TIME-TIME_min)/(25*365.25)),'--','LineWidth',3,'Color',[0.4940 0.1840 0.5560])
        plot(TIME,2.^((TIME-TIME_min)/(50*365.25)),'--','LineWidth',3,'Color','r')
        plot(TIME,2.^((TIME-TIME_min)/(10*365.25)),'--','LineWidth',3,'Color','k')
        
        han1=plot([TIME_min TIME_max],[150 150],'--k');
        
        if i==1
            ylabel('$\frac{E}{E_0}$','FontSize',30,'Interpreter','Latex');
            set(gca,'Yscale','log','ytick',[1 2 4 8 16 32 64 128 256],'YMinorTick','off');
        else
            set(gca,'Yscale','log','ytick',[1 2 4 8 16 32 64 128 256],'yticklabel',[],'YMinorTick','off');
        end
        
%         if i>1;
%             set(gca,'XTickLabel',{'' '2050' '2100' '2150' '2200' '2250' '2300'},...
%                 'xtick',datenum(2000:50:2300,1,1))
%         end
        set(gca,'XTickLabel',{'2000' '2020' '2040' '2060' '2080' '2100', ''},...
            'xtick',datenum(2000:20:2300,1,1),'FontName',FontName)
        xlabel('Year','FontName',FontName)
        
        if SCENARIO<3
            set(gca,'xticklabel',[]);
            xlabel('');
        end
        
        axis([TIME_min TIME_max 1 256]);
        
    end
    
end

annotation(gcf,'textbox',...
    [0.324750000000001 0.928420765027325 0.664624999999999 0.0633333333333344],...
    'String','Normalized increase in the probability of exceedance (i.e., "doubling")',...
    'LineStyle','none',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',20,...
    'FontName',FontName,...
    'FitBoxToText','off');

SCENARIO=1;
xp0=0.05;  yp0=0.88;
xsep1=0.225;
ysep1=0.21;

annotation(gcf,'textbox',[xp0                       yp0-(SCENARIO-1)*(height+ysep) 0.02 0.04],'String',{'A'},'FontWeight','bold','FontSize',30,'FitBoxToText','off','BackgroundColor',[0.8 0.8 0.8],'HorizontalAlignment','center','VerticalAlignment','middle','FontName','Nomada Incise');
annotation(gcf,'textbox',[xp0+xsep1+1*(width+xsep2) yp0-ysep1-(SCENARIO-1)*(height+ysep) 0.02 0.04],'String',{'B'},'FontWeight','bold','FontSize',30,'FitBoxToText','off','BackgroundColor',[0.8 0.8 0.8],'HorizontalAlignment','center','VerticalAlignment','middle','FontName','Nomada Incise');
annotation(gcf,'textbox',[xp0+xsep1+2*(width+xsep2) yp0-ysep1-(SCENARIO-1)*(height+ysep) 0.02 0.04],'String',{'C'},'FontWeight','bold','FontSize',30,'FitBoxToText','off','BackgroundColor',[0.8 0.8 0.8],'HorizontalAlignment','center','VerticalAlignment','middle','FontName','Nomada Incise');
annotation(gcf,'textbox',[xp0+xsep1+3*(width+xsep2) yp0-ysep1-(SCENARIO-1)*(height+ysep) 0.02 0.04],'String',{'D'},'FontWeight','bold','FontSize',30,'FitBoxToText','off','BackgroundColor',[0.8 0.8 0.8],'HorizontalAlignment','center','VerticalAlignment','middle','FontName','Nomada Incise');

SCENARIO=2;
annotation(gcf,'textbox',[xp0                       yp0-(SCENARIO-1)*(height+ysep) 0.02 0.04],'String',{'E'},'FontWeight','bold','FontSize',30,'FitBoxToText','off','BackgroundColor',[0.8 0.8 0.8],'HorizontalAlignment','center','VerticalAlignment','middle','FontName','Nomada Incise');
annotation(gcf,'textbox',[xp0+xsep1+1*(width+xsep2) yp0-ysep1-(SCENARIO-1)*(height+ysep) 0.02 0.04],'String',{'F'},'FontWeight','bold','FontSize',30,'FitBoxToText','off','BackgroundColor',[0.8 0.8 0.8],'HorizontalAlignment','center','VerticalAlignment','middle','FontName','Nomada Incise');
annotation(gcf,'textbox',[xp0+xsep1+2*(width+xsep2) yp0-ysep1-(SCENARIO-1)*(height+ysep) 0.02 0.04],'String',{'G'},'FontWeight','bold','FontSize',30,'FitBoxToText','off','BackgroundColor',[0.8 0.8 0.8],'HorizontalAlignment','center','VerticalAlignment','middle','FontName','Nomada Incise');
annotation(gcf,'textbox',[xp0+xsep1+3*(width+xsep2) yp0-ysep1-(SCENARIO-1)*(height+ysep) 0.02 0.04],'String',{'H'},'FontWeight','bold','FontSize',30,'FitBoxToText','off','BackgroundColor',[0.8 0.8 0.8],'HorizontalAlignment','center','VerticalAlignment','middle','FontName','Nomada Incise');

SCENARIO=3;
annotation(gcf,'textbox',[xp0                       yp0-(SCENARIO-1)*(height+ysep) 0.02 0.04],'String',{'I'},'FontWeight','bold','FontSize',30,'FitBoxToText','off','BackgroundColor',[0.8 0.8 0.8],'HorizontalAlignment','center','VerticalAlignment','middle','FontName','Nomada Incise');
annotation(gcf,'textbox',[xp0+xsep1+1*(width+xsep2) yp0-ysep1-(SCENARIO-1)*(height+ysep) 0.02 0.04],'String',{'J'},'FontWeight','bold','FontSize',30,'FitBoxToText','off','BackgroundColor',[0.8 0.8 0.8],'HorizontalAlignment','center','VerticalAlignment','middle','FontName','Nomada Incise');
annotation(gcf,'textbox',[xp0+xsep1+2*(width+xsep2) yp0-ysep1-(SCENARIO-1)*(height+ysep) 0.02 0.04],'String',{'K'},'FontWeight','bold','FontSize',30,'FitBoxToText','off','BackgroundColor',[0.8 0.8 0.8],'HorizontalAlignment','center','VerticalAlignment','middle','FontName','Nomada Incise');
annotation(gcf,'textbox',[xp0+xsep1+3*(width+xsep2) yp0-ysep1-(SCENARIO-1)*(height+ysep) 0.02 0.04],'String',{'L'},'FontWeight','bold','FontSize',30,'FitBoxToText','off','BackgroundColor',[0.8 0.8 0.8],'HorizontalAlignment','center','VerticalAlignment','middle','FontName','Nomada Incise');


% print(gcf,'figures/Figure6s_rate_of_increase_with_time_Appendix3.jpg','-djpeg','-r300');
exportgraphics(gcf,'figures/FigureA4.png','Resolution',300)


%% Exceedance Plot with Time ODDS 2

% quick check of SLR curves

% han1=figure(66); set(han1,'PaperpositionMode','auto','Position',[200 100 1600 1000],'color','w');
% 
% TEMP_PRJ=[Tmin, Tmean, Tmax]; % ????
% TIME_min=datenum(2000,1,1);
% TIME_max=datenum(2300,1,1);
% 
% for SCENARIO=1:3
%     
%     for i=1:len
%         
%         id=ID(i); % get the specific station
%         
% %         [SLR_min(i),tid]=min(S(id).SLR_med(:,SCENARIO));
%         [TEMP_min(i),tid]=min(TEMP_PRJ(:,SCENARIO));
%         
%         %plot(S(id).t_SLR,S(id).SLR_med(:,SCENARIO),'b-',S(id).t_SLR(tid),S(id).SLR_med(tid,SCENARIO),'ro'); datetick('x');
%     end
%     
%     id_neg=find(TEMP_min<-0.01);
%     
%     % S(ID(id_neg)).name
%     
%     TIME=(datenum(2000,1,1):100:datenum(2300,1,1))';
%     
%     E_TIME=NaN(length(TIME),len);
%     
%     E_TIME_DATA=NaN(length(TIME),len);
%     
%     E_TIME_GEV=NaN(length(TIME),len);
%     
%     ODDS_TIME=NaN(length(TIME),len);
%     
%     ODDS_TIME_DATA=NaN(length(TIME),len);
%     
%     ODDS_TIME_GEV=NaN(length(TIME),len);
%     
%     TR=50;
%     
% %     x1=(0:0.001:2*max(x_data))';
%     x1=(0:0.01:80)';
% %     for i=1:1:len
%     for i=1:200:len
%         
%         id=ID(i); % get the specific station
%         
%         E0=S(id).ri/TR; % get the exceedance probability of the TR-year event
%         
%         % GEV
%         mu=S(id).mu;
%         sig=S(id).sigma;
%         k=S(id).k;
%         
% %         GW=interp1(S(id).t_SLR,S(id).SLR_med(:,SCENARIO),TIME);
%         GW=interp1(ti,TEMP_PRJ(:,SCENARIO),TIME);
%         
%         %plot(S(id).t_SLR,S(id).SLR_med(:,SCENARIO),'ro',TIME,SLR,'-b.'); datetick('x');
%         
%         % empirical exceedace probability of all data
%         [Etmp,ia,~] = unique(S(id).Exceedance,'last');
%         xtmp=S(id).x(ia);
%         
%         x0=x_Tr(TR,mu,sig,k); % find the water level corresponding to the 50-year event
%         
%         Ei=interp1(xtmp,Etmp,x0,'linear');
%         
%         % empirical exceedace probability of the data of extreme events
%         len_data=length(S(id).h_max);
%         E_data=linspace(1/(len_data+1),len_data/(len_data+1),len_data)';
%         x_data=S(id).h_max;
%         
%         x0_data=interp1(E_data,x_data,E0,'linear','extrap');
%         %x0_data=interp1(E_data,x_data,E0);
%         %Ei_data=interp1(x_data,E_data,x0,'linear');
%         [tmp_x,~,~] = unique(x_data,'last');
%         tmp_y=tmp_x;
%         for jj=1:length(tmp_x)
%             ind_tmp=find(x_data==tmp_x(jj));
%             tmp_y(jj)=mean(E_data(ind_tmp));
%         end
%         
%         for kk=1:length(GW)
%             
%             % find empirical exceedance probability increase of the all values
%             E_TIME(kk,i)=interp1(xtmp+GW(kk),Etmp,x0);
%             
%             % find empirical exceedance probability increase of the extremes
% %             E_TIME_DATA(kk,i)=interp1(x_data+GW(kk),E_data,x0_data);
%             E_TIME_DATA(kk,i)=interp1(tmp_x+GW(kk),tmp_y,x0_data);
%             
%             % find GEV probability increase
%             E_TIME_GEV(kk,i)=E(x0,mu+GW(kk),sig,k);
%             
%             PLOT=0;
%             
%             if E_TIME_DATA(kk,i)<0
%                 PLOT=1;
%             end
%             
%             if PLOT
%                 
%                 figure;
%                 subplot(1,3,1);
%                 
%                 plot(xtmp,Etmp,'-b.',xtmp+GW(kk),Etmp,'-r.',...
%                     [0 max(x_data)+GW(end)],[Ei Ei],'k--',...
%                     x0,Ei,'mo',[x0 x0],[0 1],'k--',...
%                     x0,E_TIME(kk,i),'go',...
%                     [0 max(x_data)+GW(end)],[E_TIME(kk,i) E_TIME(kk,i)],'k--');
%                 
%                 axis([0 max(x_data)+GW(end) 0 1]);
%                 
%                 subplot(1,3,2);
%                 
%                 plot(x_data,E_data,'-b.',x_data+GW(kk),E_data,'-r.',...
%                     [0.9*min(x_data) max(x_data)+GW(end)],[S(id).ri/TR S(id).ri/TR],'k--',...
%                     x0_data,E0,'mo',[x0_data x0_data],[0 1],'k--',...
%                     x0_data,E_TIME_DATA(kk,i),'go',...
%                     [0.9*min(x_data) max(x_data)+GW(end)],[E_TIME_DATA(kk,i) E_TIME_DATA(kk,i)],'k--',[x0 x0],[0 1],'g--');
%                 
%                 axis([0.9*min(x_data) max(x_data)+GW(end) 0 1]);
%                 
%                 subplot(1,3,3);
%                 
%                 plot(x1,E(x1,mu,sig,k),'-b.',x1,E(x1,mu+GW(kk),sig,k),'-r.',...
%                     [0.9*min(x_data) max(x_data)+GW(end)],[S(id).ri/TR S(id).ri/TR],'k--',...
%                     x_Tr(TR,mu,sig,k),E0,'mo',[x_Tr(TR,mu,sig,k) x_Tr(TR,mu,sig,k)],[0 1],'k--',...
%                     x_Tr(TR,mu,sig,k),E_TIME_GEV(kk,i),'go',...
%                     [0.9*min(x_data) max(x_data)+GW(end)],[E_TIME_GEV(kk,i) E_TIME_GEV(kk,i)],'k--');
%                 
%                 axis([0.9*min(x_data) max(x_data)+GW(end) 0 1]);
%                 
%                 drawnow;
%                 stop
%                 
%             end
%             
%         end
%         
%         
%         
%         ODDS_TIME(:,i)=E_TIME(:,i)./(1-E_TIME(:,i));
%         ODDS_TIME_DATA(:,i)=E_TIME_DATA(:,i)./(1-E_TIME_DATA(:,i));
%         ODDS_TIME_GEV(:,i)=E_TIME_GEV(:,i)./(1-E_TIME_GEV(:,i));
%         
%         ODDSi=Ei./(1-Ei);
%         ODDS0=E0./(1-E0);
%         
%         color1=cmap(idx(i),:); % color based on clusters
%         
%         if ismember(i,id_neg)
%             color1='k';           % single color
%         end
%         
%         ALPHA=0.15;
%         
%         xp0=0.045;
%         yp0=0.66;
%         width=0.205;
%         height=0.27;
%         xsep1=0.05;
%         xsep2=0.02;
%         ysep=0.03;
%         
% 
%         
%         subplot('position',[xp0+xsep1+1*(width+xsep2) yp0-(SCENARIO-1)*(height+ysep) width height]); hold on; box on;
%         han1=plot(TIME,ODDS_TIME(:,i)/ODDSi,'Color',color1); han1.Color(4)=ALPHA;  % set transparency
%         set(han1,'LineWidth',2);
%         if SCENARIO==1
%             title({'for all data'});
%         end
%         set(gca,'FontSize',14);
%         datetick('x');
% 
%         subplot('position',[xp0+xsep1+2*(width+xsep2) yp0-(SCENARIO-1)*(height+ysep) width height]); hold on; box on;
%         han1=plot(TIME,ODDS_TIME_DATA(:,i)/ODDS0,'Color',color1); han1.Color(4)=ALPHA;  % set transparency
%         set(han1,'LineWidth',2);
%         if SCENARIO==1;
%             title({'for extreme events data'});
%         end
%         set(gca,'FontSize',14);
%         datetick('x');
%         
%         subplot('position',[xp0+xsep1+3*(width+xsep2) yp0-(SCENARIO-1)*(height+ysep) width height]); hold on; box on;
%         han1=plot(TIME,ODDS_TIME_GEV(:,i)/ODDS0,'Color',color1); han1.Color(4)=ALPHA;  % set transparency
%         set(han1,'LineWidth',2);
%         if SCENARIO==1;
%             title({'for GEV model of extremes'})
%         end
%         set(gca,'FontSize',14);
%         datetick('x');
%         
%     end
%     
%     
%     subplot('position',[xp0 yp0-(SCENARIO-1)*(height+ysep) width height]); hold on; box on;
%     han1=plot(TIME,GW,'Color','k'); %han1.Color(4)=ALPHA;  % set transparency
%     set(han1,'LineWidth',2);
%     if SCENARIO==1
%         title({'ensemble mean GTR'});
%         ylabel('Best temp. scenario [C]','FontSize',14);
%     elseif SCENARIO==2
%         ylabel('Average temp. scenario [C]','FontSize',14);
%     elseif SCENARIO==3
%         ylabel('Worst temp. scenario [C]','FontSize',14);
%     end
%     
%     axis([TIME_min TIME_max -0.5 15]);
%     
%     
%     set(gca,'FontSize',14);
%     datetick('x');
%     if SCENARIO<3
%         set(gca,'xticklabel',[]);
%     end
%     
%     
%     TIME_min=datenum(2000,1,1);
%     TIME_max=datenum(2300,1,1);
%     
%     for i=1:3
%         
%         subplot('position',[xp0+xsep1+i*(width+xsep2) yp0-(SCENARIO-1)*(height+ysep) width height]); hold on; box on;
% %         
% %         text(731307.357995124,75.4065102899739,'5 yr','Rotation',85,'FontSize',14,'color','r');
% %         text(739816.264607017,66.8559432854146 ,'10 yrs','Rotation',65,'FontSize',14,'color','k');
% %         text(760204.522935306,12.2468383865414,'25 yrs','Rotation',24,'FontSize',14,'color','g','BackgroundColor',[1 1 1]);
%         
%         text(737986.09325084, 73.303021521361,'5 yr','Rotation',85,'FontSize',14,'color','r');
%         text(746077.578529558, 48.9760389940741 ,'10 yrs','Rotation',70,'FontSize',14,'color','k');
%         text(773144.572316924, 16.2514969684399,'25 yrs','Rotation',50,'FontSize',14,'color','b');%,'BackgroundColor',[1 1 1]
%         
%         han1=plot(TIME,2.^((TIME-TIME_min)/(10*365.25)),'--k',TIME,2.^((TIME-TIME_min)/(365.25*5)),'--r',TIME,2.^((TIME-TIME_min)/(25*365.25)),'--b');
%         set(han1,'LineWidth',3);
%         
% %         han1=plot([TIME_min TIME_max],[150 150],'--k');
%         
%         if i==1;
% %             ylabel('$\frac{Odds}{Odds_0}$','FontSize',30,'Interpreter','Latex');
%             ylabel('$\frac{\mathrm{Odds}}{\mathrm{Odds}_0}$','FontSize',30,'Interpreter','Latex');
%             set(gca,'Yscale','log','ytick',[1 2 4 8 16 32 64 128 256],'YMinorTick','off');
%         else
%             set(gca,'Yscale','log','ytick',[1 2 4 8 16 32 64 128 256],'yticklabel',[],'YMinorTick','off');
%         end
%         
%         if i>1;
%             set(gca,'XTickLabel',{'' '2050' '2100' '2150' '2200' '2250' '2300'},...
%                 'xtick',datenum(2000:50:2300,1,1))
%             
%         end
%         
%         if SCENARIO<3
%             set(gca,'xticklabel',[]);
%         end
%         
%         axis([TIME_min TIME_max 1 256]);
%         
%     end
%     
% end
% 
% annotation(gcf,'textbox',[0.454750000000001 0.938666666666667 0.607333333333333 0.0633333333333344],'String',{'Normalized odds increase (i.e., "doubling")'},'FitBoxToText','off','linestyle','none','FontWeight','bold','FontSize',16);
% 
% % annotation(gcf,'textbox',[0.340500000000001 0.825000000000001 0.0236666666666666 0.0583333333333342],'String',{'A'},'FontWeight','bold','FontSize',16,'FitBoxToText','off','BackgroundColor',[1 1 1]);
% % annotation(gcf,'textbox',[0.340500000000001+1*(width+xsep) 0.825000000000001 0.0236666666666666 0.0583333333333342],'String',{'B'},'FontWeight','bold','FontSize',16,'FitBoxToText','off','BackgroundColor',[1 1 1]);
% % annotation(gcf,'textbox',[0.340500000000001+2*(width+xsep) 0.825000000000001 0.0236666666666666 0.0583333333333342],'String',{'C'},'FontWeight','bold','FontSize',16,'FitBoxToText','off','BackgroundColor',[1 1 1]);
% 
% print(gcf,'figures/Figure66s_rate_of_increase_with_time.jpg','-djpeg','-r300');


%% Exceedance Plot with Time ODDS - Required for fig 7 (Appendix 4)

lat=NaN(len,1);
lon=NaN(len,1);

clus=NaN(len,1);

tau_2x=NaN(len,3);
tau_2x_DATA=NaN(len,3);
tau_2x_GEV=NaN(len,3);

tau_2x_2050=NaN(len,3);
tau_2x_2050_DATA=NaN(len,3);
tau_2x_2050_GEV=NaN(len,3);

tau_2x_fit=NaN(len,3);
tau_2x_fit_DATA=NaN(len,3);
tau_2x_fit_GEV=NaN(len,3);

TIME_min=datenum(2000,1,1);
% TIME_max=datenum(2300,1,1);
TIME_max=datenum(2100,1,1);
TEMP_PRJ=[Tmin, Tmean, Tmax]; % Low-end, Mean, High-end RCP 8.5 

% quick check of SLR curves

han1=figure(66); set(han1,'PaperpositionMode','auto','Position',[200 100 1600 1000],'color','w');

for SCENARIO=1:3 
    
    for i=1:len
        
%         id=ID(i); % get the specific station
%         [SLR_min(i),tid]=min(S(id).SLR_med(:,SCENARIO));
        [TEMP_min(i),tid]=min(TEMP_PRJ(:,SCENARIO));
        
        %plot(S(id).t_SLR,S(id).SLR_med(:,SCENARIO),'b-',S(id).t_SLR(tid),S(id).SLR_med(tid,SCENARIO),'ro'); datetick('x');
    end
    
    id_neg=find(TEMP_min<-0.01);
    
    % S(ID(id_neg)).name
    
    
%     TIME=(TIME0:100:datenum(2100,1,1))';
%   New changes: 
%     TIME0=datenum(2000,1,1); % Global
%     TIME=(TIME0:100:datenum(2150,1,1))'; % Global
    TIME0=T_proj_local.time_start; % Local
    TIME=(TIME0:100:datenum(2100,1,1))'; % Global
    
    TIME_threshold1=datenum(2000,1,1);        % SLR threshold for fit of tau 2x
    TIME_threshold2=datenum(2050,1,1);        % SLR threshold for fit of tau 2x
    
    idsel1=TIME>=TIME_threshold1 & TIME<=TIME_threshold2;
    
    TIME_threshold1=datenum(2025,1,1);        % SLR threshold for fit of tau 2x
    TIME_threshold2=datenum(2075,1,1);        % SLR threshold for fit of tau 2x
    
    idsel2=TIME>=TIME_threshold1 & TIME<=TIME_threshold2;
    
    E_TIME=NaN(length(TIME),len);
    
    E_TIME_DATA=NaN(length(TIME),len);
    
    E_TIME_GEV=NaN(length(TIME),len);
    
    ODDS_TIME=NaN(length(TIME),len);
    ODDS_TIME_DATA=NaN(length(TIME),len);
    ODDS_TIME_GEV=NaN(length(TIME),len);
    
    TR=50;
    
    x1=(0:0.001:2*max(x_data))';
    
%     for i=1:len
    for i=1:1:len % This will not go over all the stations
        [SCENARIO i]
        
        id=ID(i); % get the specific station
        
        E0=S(id).ri/TR; % get the exceedance probability of the TR-year event
        
        % GEV
        mu=S(id).mu;
        sig=S(id).sigma;
        k=S(id).k;
        
        lat(i)=S(id).lat;
        lon(i)=S(id).lon;
        
        % cluster
        clus(i)=idx(i);
         
%         SLR=interp1(S(id).t_SLR,S(id).SLR_med(:,SCENARIO),TIME);
%       New change:
%         SLR=interp1(ti, TEMP_PRJ(:,SCENARIO),TIME); % Global
        SLR=interp1(S(id).ti, S(id).TR_scenario(:,SCENARIO),TIME); % Local
        
        %plot(S(id).t_SLR,S(id).SLR_med(:,SCENARIO),'ro',TIME,SLR,'-b.'); datetick('x');
        
        % empirical exceedace probability of all data
        [Etmp,ia,~] = unique(S(id).Exceedance,'last');
        xtmp=S(id).x(ia);
        
        x0=x_Tr(TR,mu,sig,k); % find the water level corresponding to the 50-year event
        
        Ei=interp1(xtmp,Etmp,x0,'linear');
        
        % empirical exceedace probability of the data of extreme events
        len_data=length(S(id).h_max);
        E_data=linspace(1/(len_data+1),len_data/(len_data+1),len_data)';
        x_data=S(id).h_max;
        
        x0_data=interp1(E_data,x_data,E0,'linear','extrap');
        %Ei_data=interp1(x_data,E_data,x0,'linear');
        
        for kk=1:length(SLR)
            
            % find empirical exceedance probability increase of the all values
            E_TIME(kk,i)=interp1(xtmp+SLR(kk),Etmp,x0);
            
            % find empirical exceedance probability increase of the extremes
            [~,ind]=unique(x_data); % ???? I added. Is this Okay?
            E_TIME_DATA(kk,i)=interp1(x_data(ind)+SLR(kk),E_data(ind),x0_data);
            
            % find GEV probability increase
            E_TIME_GEV(kk,i)=E(x0,mu+SLR(kk),sig,k);
            
            PLOT=0;
            
            if E_TIME_DATA(kk,i)<0
                PLOT=1;
            end
            
            if PLOT
                
                figure;
                subplot(1,3,1);
                
                plot(xtmp,Etmp,'-b.',xtmp+SLR(kk),Etmp,'-r.',...
                    [0 max(x_data)+SLR(end)],[Ei Ei],'k--',...
                    x0,Ei,'mo',[x0 x0],[0 1],'k--',...
                    x0,E_TIME(kk,i),'go',...
                    [0 max(x_data)+SLR(end)],[E_TIME(kk,i) E_TIME(kk,i)],'k--');
                
                axis([0 max(x_data)+SLR(end) 0 1]);
                
                subplot(1,3,2);
                
                plot(x_data,E_data,'-b.',x_data+SLR(kk),E_data,'-r.',...
                    [0.9*min(x_data) max(x_data)+SLR(end)],[S(id).ri/TR S(id).ri/TR],'k--',...
                    x0_data,E0,'mo',[x0_data x0_data],[0 1],'k--',...
                    x0_data,E_TIME_DATA(kk,i),'go',...
                    [0.9*min(x_data) max(x_data)+SLR(end)],[E_TIME_DATA(kk,i) E_TIME_DATA(kk,i)],'k--',[x0 x0],[0 1],'g--');
                
                axis([0.9*min(x_data) max(x_data)+SLR(end) 0 1]);
                
                subplot(1,3,3);
                
                plot(x1,E(x1,mu,sig,k),'-b.',x1,E(x1,mu+SLR(kk),sig,k),'-r.',...
                    [0.9*min(x_data) max(x_data)+SLR(end)],[S(id).ri/TR S(id).ri/TR],'k--',...
                    x_Tr(TR,mu,sig,k),E0,'mo',[x_Tr(TR,mu,sig,k) x_Tr(TR,mu,sig,k)],[0 1],'k--',...
                    x_Tr(TR,mu,sig,k),E_TIME_GEV(kk,i),'go',...
                    [0.9*min(x_data) max(x_data)+SLR(end)],[E_TIME_GEV(kk,i) E_TIME_GEV(kk,i)],'k--');
                
                axis([0.9*min(x_data) max(x_data)+SLR(end) 0 1]);
                
                drawnow;
                stop
                
            end
            
        end
        
        ODDS_TIME(:,i)=E_TIME(:,i)./(1-E_TIME(:,i));
        ODDS_TIME_DATA(:,i)=E_TIME_DATA(:,i)./(1-E_TIME_DATA(:,i));
        ODDS_TIME_GEV(:,i)=E_TIME_GEV(:,i)./(1-E_TIME_GEV(:,i));
        
        ODDSi=Ei./(1-Ei);
        ODDS0=E0./(1-E0);
        
        color1=cmap(idx(i),:); % color based on clusters
        
        if ismember(i,id_neg)
            color1='k';           % single color
        end
        
        ALPHA=0.15;
        
        xp0=0.045;
        yp0=0.66;
        width=0.205;
        height=0.27;
        xsep1=0.05;
        xsep2=0.025;
        ysep=0.03;
        
        
        subplot('position',[xp0 yp0-(SCENARIO-1)*(height+ysep) width height]); hold on; box on;
        han1=plot(TIME,SLR,'Color', color1); han1.Color(4)=ALPHA;  % set transparency
        set(han1,'LineWidth',2);
        
        
        subplot('position',[xp0+xsep1+1*(width+xsep2) yp0-(SCENARIO-1)*(height+ysep) width height]); hold on; box on;
        han1=plot(TIME,ODDS_TIME(:,i)/ODDSi,'Color',color1); han1.Color(4)=ALPHA;  % set transparency
        set(han1,'LineWidth',2);

        
        datetick('x');
        if SCENARIO==1
            title({'All data'},'FontSize',18,'FontName',FontName);
        end
        set(gca,'FontSize',14,'FontName',FontName);

        subplot('position',[xp0+xsep1+2*(width+xsep2) yp0-(SCENARIO-1)*(height+ysep) width height]); hold on; box on;
        han1=plot(TIME,ODDS_TIME_DATA(:,i)/ODDS0,'Color',color1); han1.Color(4)=ALPHA;  % set transparency
        set(han1,'LineWidth',2);
   
        
        datetick('x');
        if SCENARIO==1
            title({'Extreme events data'},'FontSize',18,'FontName',FontName);
        end
        set(gca,'FontSize',14,'FontName',FontName);
        
        subplot('position',[xp0+xsep1+3*(width+xsep2) yp0-(SCENARIO-1)*(height+ysep) width height]); hold on; box on;
        han1=plot(TIME,ODDS_TIME_GEV(:,i)/ODDS0,'Color',color1); han1.Color(4)=ALPHA;  % set transparency
        set(han1,'LineWidth',2);

        
        datetick('x');
        if SCENARIO==1
            title({'GEV model from extreme events data'},'FontSize',18,'FontName',FontName);
        end
        set(gca,'FontSize',14,'FontName',FontName);
        
        % find the time scale to double the odds

        xtmp=TIME(idsel1);
        
        % find the amount of SLR to double odds
        
        ytmp=log2(ODDS_TIME(idsel1,i)/ODDSi);      idn=~isnan(ytmp); p=polyfit(xtmp(idn),ytmp(idn),1);   % fit straight line
        
        tau_2x(i,SCENARIO)=1/(nanmean(diff(ytmp)./diff(xtmp)));                 % find tau 2x
        tau_2x_fit(i,SCENARIO)=1/p(1);                                          % find tau 2x using a regression
        
        ytmp=log2(ODDS_TIME_DATA(idsel1,i)/ODDSi);      idn=~isnan(ytmp); p=polyfit(xtmp(idn),ytmp(idn),1);   % fit straight line
        
        tau_2x_DATA(i,SCENARIO)=1/(nanmean(diff(ytmp)./diff(xtmp)));            % find tau 2x
        tau_2x_fit_DATA(i,SCENARIO)=1/p(1);                                     % find tau 2x using a regression
        
        ytmp=log2(ODDS_TIME_GEV(idsel1,i)/ODDSi);      idn=~isnan(ytmp); p=polyfit(xtmp(idn),ytmp(idn),1);   % fit straight line
        
        tau_2x_GEV(i,SCENARIO)=1/(nanmean(diff(ytmp)./diff(xtmp)));             % find tau 2x
        tau_2x_fit_GEV(i,SCENARIO)=1/p(1);                                      % find tau 2x using a regression
        
        % second time frame
        
        xtmp=TIME(idsel2);
        
        % find the amount of SLR to double odds
        
        ytmp=log2(ODDS_TIME(idsel2,i)/ODDSi);  % fit straight line
        
        tau_2x_2050(i,SCENARIO)=1/(nanmean(diff(ytmp)./diff(xtmp)));                 % find tau 2x
        
        ytmp=log2(ODDS_TIME_DATA(idsel2,i)/ODDSi);   % fit straight line
        
        tau_2x_2050_DATA(i,SCENARIO)=1/(nanmean(diff(ytmp)./diff(xtmp)));            % find tau 2x
        
        ytmp=log2(ODDS_TIME_GEV(idsel2,i)/ODDSi);      % fit straight line
        
        tau_2x_2050_GEV(i,SCENARIO)=1/(nanmean(diff(ytmp)./diff(xtmp)));             % find tau 2x
       
        if isnan(ODDSi)          % if the initial odds are undefined because the emperical probabiliy of the 50-year event cannot be interpolated... set sigma 2x to be a NaN
            tau_2x(i,SCENARIO)=NaN;
            tau_2x_fit(i,SCENARIO)=NaN;
            tau_2x_2050(i,SCENARIO)=NaN;
        end
        
        PLOT3=0;
        if PLOT3
            
            figure(667); hold on; box on;
            
            subplot(1,3,1);
            
            han1=plot(TIME,ODDS_TIME(:,i)/ODDSi,'b',TIME,2.^p1(2)*2.^(TIME/tau_2x(i,SCENARIO)),'r');
            set(han1,'LineWidth',2);
            set(gca,'Yscale','log','ytick',[1 2 4 8 16 32 64 128 256 512 1024],'YMinorTick','off');  datetick('x')
            axis([TIME_min TIME_max 1 1024]);
            
            subplot(1,3,2);
            
            han1=plot(TIME,ODDS_TIME_DATA(:,i)/ODDS0,'b',TIME,2.^p2(2)*2.^(TIME/tau_2x_DATA(i,SCENARIO)),'r');
            set(han1,'LineWidth',2);
            set(gca,'Yscale','log','ytick',[1 2 4 8 16 32 64 128 256 512 1024],'YMinorTick','off');  datetick('x')
            axis([TIME_min TIME_max 1 1024]);
            
            subplot(1,3,3);
            
            han1=plot(TIME,ODDS_TIME_GEV(:,i)/ODDS0,'b',TIME,2.^p3(2)*2.^(TIME/tau_2x_GEV(i,SCENARIO)),'r');
            set(han1,'LineWidth',2);
            set(gca,'Yscale','log','ytick',[1 2 4 8 16 32 64 128 256 512 1024],'YMinorTick','off'); datetick('x')
            axis([TIME_min TIME_max 1 1024]);
           
            drawnow;
            
        end
        
    end
    
    
    subplot('position',[xp0 yp0-(SCENARIO-1)*(height+ysep) width height]); hold on; box on;
%     han1=plot(TIME,SLR,'Color', 'b'); %han1.Color(4)=ALPHA;  % set transparency
%     set(han1,'LineWidth',2);
%     if SCENARIO==1
%         title({'ensemble mean GTR'});
%         ylabel('Best temp. scenario [C]','FontSize',14);
%     elseif SCENARIO==2
%         ylabel('Average temp. scenario [C]','FontSize',14);
%     elseif SCENARIO==3
%         ylabel('Worst temp. scenario [C]','FontSize',14);
%     end

    if SCENARIO==1
        title({'Local temperature projections'},'FontSize',16,'FontName',FontName);
        ylabel('Low-end RCP 8.5 TR [\circC]','FontSize',14);
    elseif SCENARIO==2
        ylabel('Intermediate RCP 8.5 TR [\circC]','FontSize',14);
    elseif SCENARIO==3
        ylabel('High-end RCP 8.5 TR [\circC]','FontSize',14);
    end

%     axis([TIME_min TIME_max -0.5 15]);
    axis([TIME_min TIME_max -0.5 10.5]);

    set(gca,'FontSize',14,'FontName',FontName);
    datetick('x');
    xlabel('Year');
    if SCENARIO<3   
        set(gca,'xticklabel',[]);
        xlabel('');
    end

    if i == i(end)
        if SCENARIO == 1
            AVG_Global = Tmin(1:1212,:);
            AVG_Local = mean(T_low_stations(:,1:101),1);
        elseif SCENARIO == 2
            AVG_Global = Tmean(1:1212,:);
            AVG_Local = mean(T_mean_stations(:,1:101),1);
        elseif SCENARIO == 3
            AVG_Global = Tmax(1:1212,:);
            AVG_Local = mean(T_high_stations(:,1:101),1);
        end

        AVG_Global_int = interp1(1:length(AVG_Global),AVG_Global',1:length(AVG_Global)/length(TIME):length(AVG_Global));
        AVG_Local_int = interp1(1:length(AVG_Local),AVG_Local,linspace(1,length(AVG_Local),length(TIME)));

        han11 = plot(TIME,AVG_Global_int,'--k','LineWidth',2);
        han22 = plot(TIME,AVG_Local_int,'--b','LineWidth',2);

        if SCENARIO == 1
            legend([han11 han22],'Global average','North America average','FontSize',14,'FontName',FontName);
        end
    end


    TIME_min=datenum(2000,1,1);
    TIME_max=datenum(2100,1,1);
    
    for i=1:3
        
        subplot('position',[xp0+xsep1+i*(width+xsep2) yp0-(SCENARIO-1)*(height+ysep) width height]); hold on; box on;
        
%         text(731307.357995124,75.4065102899739,'1 yr','Rotation',85,'FontSize',14,'color','r');
%         text(739816.264607017,66.8559432854146 ,'5 yrs','Rotation',65,'FontSize',14,'color','k');
%         text(760204.522935306,12.2468383865414,'25 yrs','Rotation',24,'FontSize',14,'color','g','BackgroundColor',[1 1 1]);

        text(731306.861571456, 82.6691307776895,'1 year','Rotation',88,'FontSize',14,'color','m','FontName',FontName);
        text(739435.75474142, 59.3537289020805 ,'5 years','Rotation',67,'FontSize',14,'color','b','FontName',FontName);
        text(746790.844100608, 33.7572473193161 ,'10 years','Rotation',50,'FontSize',14,'color','k','FontName',FontName);
        text(755415.747120308, 8.80139136611196 ,'25 years','Rotation',25,'FontSize',14,'color',[0.4940 0.1840 0.5560],'FontName',FontName);
        text(757698.809620308, 3.64269028542452 ,'50 years','Rotation',15,'FontSize',14,'color','r','FontName',FontName);

%         han1=plot(TIME,2.^((TIME-TIME_min)/(5*365.25)),'--k',TIME,2.^((TIME-TIME_min)/365.25),'--r',TIME,2.^((TIME-TIME_min)/(25*365.25)),'--g');

        plot(TIME,2.^((TIME-TIME_min)/(5*365.25)),'--','LineWidth',3,'Color','b'); hold on
        plot(TIME,2.^((TIME-TIME_min)/365.25),'--','LineWidth',3,'Color','m')
        plot(TIME,2.^((TIME-TIME_min)/(25*365.25)),'--','LineWidth',3,'Color',[0.4940 0.1840 0.5560])
        plot(TIME,2.^((TIME-TIME_min)/(50*365.25)),'--','LineWidth',3,'Color','r')
        plot(TIME,2.^((TIME-TIME_min)/(10*365.25)),'--','LineWidth',3,'Color','k')

%         han1=plot([TIME_min TIME_max],[150 150],'--k');
       

        
%         if i>1
%             set(gca,'XTickLabel',{'' '2020' '2040' '2060' '2080' '2100'})
%         end
        set(gca,'XTickLabel',{'2000' '2020' '2040' '2060' '2080' '2100', ''},...
            'xtick',datenum(2000:20:2300,1,1),'FontSize',14,'FontName',FontName)
        xlabel('Year');

        if SCENARIO<3   
            set(gca,'xticklabel',[]);
            xlabel('');
        end

        if i==1
            ylabel('$\frac{\mathrm{Odds}}{\mathrm{Odds}_0}$','FontSize',30,'Interpreter','latex');
            set(gca,'Yscale','log','ytick',[1 2 4 8 16 32 64 128 256],'YMinorTick','off');
        else
            set(gca,'Yscale','log','ytick',[1 2 4 8 16 32 64 128 256],'yticklabel',[],'YMinorTick','off');
        end
        
        axis([TIME_min TIME_max 1 256]);
        
    end
    
end

annotation(gcf,'textbox',...
    [0.324750000000001 0.928420765027325 0.664624999999999 0.0633333333333344],...
    'String','Normalized increase in the odds of exceedance (i.e., "doubling")',...
    'LineStyle','none',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',20,...
    'FontName',FontName,...
    'FitBoxToText','off');

SCENARIO=1;
xp0=0.05;  yp0=0.88;
xsep1=0.225;
ysep1=0.21;

annotation(gcf,'textbox',[xp0                       yp0-(SCENARIO-1)*(height+ysep) 0.02 0.04],'String',{'A'},'FontWeight','bold','FontSize',30,'FitBoxToText','off','BackgroundColor',[0.8 0.8 0.8],'HorizontalAlignment','center','VerticalAlignment','middle','FontName','Nomada Incise');
annotation(gcf,'textbox',[xp0+xsep1+1*(width+xsep2) yp0-ysep1-(SCENARIO-1)*(height+ysep) 0.02 0.04],'String',{'B'},'FontWeight','bold','FontSize',30,'FitBoxToText','off','BackgroundColor',[0.8 0.8 0.8],'HorizontalAlignment','center','VerticalAlignment','middle','FontName','Nomada Incise');
annotation(gcf,'textbox',[xp0+xsep1+2*(width+xsep2) yp0-ysep1-(SCENARIO-1)*(height+ysep) 0.02 0.04],'String',{'C'},'FontWeight','bold','FontSize',30,'FitBoxToText','off','BackgroundColor',[0.8 0.8 0.8],'HorizontalAlignment','center','VerticalAlignment','middle','FontName','Nomada Incise');
annotation(gcf,'textbox',[xp0+xsep1+3*(width+xsep2) yp0-ysep1-(SCENARIO-1)*(height+ysep) 0.02 0.04],'String',{'D'},'FontWeight','bold','FontSize',30,'FitBoxToText','off','BackgroundColor',[0.8 0.8 0.8],'HorizontalAlignment','center','VerticalAlignment','middle','FontName','Nomada Incise');

SCENARIO=2;
annotation(gcf,'textbox',[xp0                       yp0-(SCENARIO-1)*(height+ysep) 0.02 0.04],'String',{'E'},'FontWeight','bold','FontSize',30,'FitBoxToText','off','BackgroundColor',[0.8 0.8 0.8],'HorizontalAlignment','center','VerticalAlignment','middle','FontName','Nomada Incise');
annotation(gcf,'textbox',[xp0+xsep1+1*(width+xsep2) yp0-ysep1-(SCENARIO-1)*(height+ysep) 0.02 0.04],'String',{'F'},'FontWeight','bold','FontSize',30,'FitBoxToText','off','BackgroundColor',[0.8 0.8 0.8],'HorizontalAlignment','center','VerticalAlignment','middle','FontName','Nomada Incise');
annotation(gcf,'textbox',[xp0+xsep1+2*(width+xsep2) yp0-ysep1-(SCENARIO-1)*(height+ysep) 0.02 0.04],'String',{'G'},'FontWeight','bold','FontSize',30,'FitBoxToText','off','BackgroundColor',[0.8 0.8 0.8],'HorizontalAlignment','center','VerticalAlignment','middle','FontName','Nomada Incise');
annotation(gcf,'textbox',[xp0+xsep1+3*(width+xsep2) yp0-ysep1-(SCENARIO-1)*(height+ysep) 0.02 0.04],'String',{'H'},'FontWeight','bold','FontSize',30,'FitBoxToText','off','BackgroundColor',[0.8 0.8 0.8],'HorizontalAlignment','center','VerticalAlignment','middle','FontName','Nomada Incise');

SCENARIO=3;
annotation(gcf,'textbox',[xp0                       yp0-(SCENARIO-1)*(height+ysep) 0.02 0.04],'String',{'I'},'FontWeight','bold','FontSize',30,'FitBoxToText','off','BackgroundColor',[0.8 0.8 0.8],'HorizontalAlignment','center','VerticalAlignment','middle','FontName','Nomada Incise');
annotation(gcf,'textbox',[xp0+xsep1+1*(width+xsep2) yp0-ysep1-(SCENARIO-1)*(height+ysep) 0.02 0.04],'String',{'J'},'FontWeight','bold','FontSize',30,'FitBoxToText','off','BackgroundColor',[0.8 0.8 0.8],'HorizontalAlignment','center','VerticalAlignment','middle','FontName','Nomada Incise');
annotation(gcf,'textbox',[xp0+xsep1+2*(width+xsep2) yp0-ysep1-(SCENARIO-1)*(height+ysep) 0.02 0.04],'String',{'K'},'FontWeight','bold','FontSize',30,'FitBoxToText','off','BackgroundColor',[0.8 0.8 0.8],'HorizontalAlignment','center','VerticalAlignment','middle','FontName','Nomada Incise');
annotation(gcf,'textbox',[xp0+xsep1+3*(width+xsep2) yp0-ysep1-(SCENARIO-1)*(height+ysep) 0.02 0.04],'String',{'L'},'FontWeight','bold','FontSize',30,'FitBoxToText','off','BackgroundColor',[0.8 0.8 0.8],'HorizontalAlignment','center','VerticalAlignment','middle','FontName','Nomada Incise');


% print(gcf,'figures/Figure6s_rate_of_increase_with_time_ODDS_Appendix4.jpg','-djpeg','-r300');
exportgraphics(gcf,'figures/FigureA5.png','Resolution',300)
% notdoneyet

tau_2x=tau_2x/365.25;                  % find tau 2x
tau_2x_DATA=tau_2x_DATA/365.25;
tau_2x_GEV=tau_2x_GEV/365.25;

tau_2x_fit=tau_2x_fit/365.25;                  % find tau 2x
tau_2x_fit_DATA=tau_2x_fit_DATA/365.25;
tau_2x_fit_GEV=tau_2x_fit_GEV/365.25;

tau_2x_2050=tau_2x_2050/365.25;                  % find tau 2x
tau_2x_2050_DATA=tau_2x_2050_DATA/365.25;
tau_2x_2050_GEV=tau_2x_2050_GEV/365.25;

tau_2x(tau_2x<=0)=NaN;
tau_2x_DATA(tau_2x_DATA<=0)=NaN;    
tau_2x_GEV(tau_2x_GEV<=0)=NaN;  

tau_2x(tau_2x>100)=NaN;
tau_2x_DATA(tau_2x_DATA>100)=NaN;    
tau_2x_GEV(tau_2x_GEV>100)=NaN;  

tau_2x_fit(tau_2x<=0)=NaN;
tau_2x_fit_DATA(tau_2x_DATA<=0)=NaN;    
tau_2x_fit_GEV(tau_2x_GEV<=0)=NaN;  

tau_2x_fit(tau_2x>100)=NaN;
tau_2x_fit_DATA(tau_2x_DATA>100)=NaN;    
tau_2x_fit_GEV(tau_2x_GEV>100)=NaN;  

tau_2x_2050(tau_2x<=0)=NaN;
tau_2x_2050_DATA(tau_2x_DATA<=0)=NaN;    
tau_2x_2050_GEV(tau_2x_GEV<=0)=NaN;  

tau_2x_2050(tau_2x>100)=NaN;
tau_2x_2050_DATA(tau_2x_DATA>100)=NaN;    
tau_2x_2050_GEV(tau_2x_GEV>100)=NaN;  

clear SLR



%% plot tau 2x (Figure 7)

SCENARIO=2; % Mean RCP 8.5

xp1=0.02;
yp1=0.08;
width=0.6; width2=0.32;
height=0.42;
height2=0.39;
xsep=0.03;
ysep=0.03;

han1=figure(32); set(han1,'PaperpositionMode','auto','Position',[200 100 1000 700],'color','w'); hold on; box on;

tau=tau_2x;
id=~isnan(tau(:,SCENARIO));
%Subplot A
subplot('position',[xp1 yp1+height+ysep width height]); hold on; box on;
geoshow('landareas.shp', 'FaceColor',[0.6 0.6 0.6],'EdgeColor',[0.6 0.6 0.6]); 
axis equal; axis tight;
scatter(lon(id),lat(id),30*ones(size(lat(id))),tau(id,SCENARIO),'filled'); colormap(gca,flipud(jet)); caxis([5 35]); colorbar;
axis([-180 -60 0 80]);

set(gca,'FontSize',14,'XTick',[-180:20:-60],'YTick',[0:10:80],'Layer','top','FontName',FontName);
ax = gca;
% Simply color an XTickLabel
ax.XTickLabel{1} = ['\color{white}' ax.XTickLabel{1}];
ax.XTickLabel{2} = ['\color{white}' ax.XTickLabel{2}];
ax.XTickLabel{3} = ['\color{white}' ax.XTickLabel{3}];
ax.XTickLabel{4} = ['\color{white}' ax.XTickLabel{4}];
ax.XTickLabel{5} = ['\color{white}' ax.XTickLabel{5}];
ax.XTickLabel{6} = ['\color{white}' ax.XTickLabel{6}];
ax.XTickLabel{7} = ['\color{white}' ax.XTickLabel{7}];

ylabel('2000-2050','FontSize',28)
title('Spatial variation in \tau','FontSize',18,'FontName',FontName);

tau=tau_2x_2050;
id=~isnan(tau(:,SCENARIO));
%Subplot C
subplot('position',[xp1 yp1 width height]); hold on; box on;
geoshow('landareas.shp', 'FaceColor',[0.6 0.6 0.6],'EdgeColor',[0.6 0.6 0.6]); 
axis equal; axis tight;
scatter(lon(id),lat(id),30*ones(size(lat(id))),tau(id,SCENARIO),'filled'); colormap(gca,flipud(jet)); caxis([5 35]); colorbar;
axis([-180 -60 0 80]);
set(gca,'FontSize',14,'XTick',[-180:20:-60],'YTick',[0:10:80],'Layer','top','FontName',FontName);
ylabel('2025-2075','FontSize',28)
tau=tau_2x;

tau_max=26;
count_max=70;
%Subplot B
subplot('position',[xp1+width+xsep yp1+height+ysep width2 height]); hold on; box on;
xe=(0:2:52)';

ybar=NaN(length(xe),4);

cmap2=cmap([2 3 1 4],:);

for i=1:length(xe)-1
   ybar(i,1)=sum(idx==2 & tau(:,SCENARIO)>=xe(i) & tau(:,SCENARIO)<xe(i+1));
   ybar(i,2)=sum(idx==3 & tau(:,SCENARIO)>=xe(i) & tau(:,SCENARIO)<xe(i+1));
   ybar(i,3)=sum(idx==1 & tau(:,SCENARIO)>=xe(i) & tau(:,SCENARIO)<xe(i+1));
   ybar(i,4)=sum(idx==4 & tau(:,SCENARIO)>=xe(i) & tau(:,SCENARIO)<xe(i+1));
end

ybar(end,1)=sum(idx==2 & tau(:,SCENARIO)>xe(end));
ybar(end,2)=sum(idx==3 & tau(:,SCENARIO)>xe(end));
ybar(end,3)=sum(idx==1 & tau(:,SCENARIO)>xe(end));
ybar(end,4)=sum(idx==4 & tau(:,SCENARIO)>xe(end));


bh=bar(xe,ybar,'stacked'); 
% colormap(gca,cmap2);  
for i=1:4
    bh(i).FaceColor = cmap2(i,:);
end


ylabel('Count (No. of stations)','FontSize',14);
set(gca,'FontSize',14,...
    'XLim',[5 35],'XTick',[0:5:35],...
    'YLim',[0 1250],...
    'FontName',FontName,...
    'Layer','top');
% 'YLim',[0 count_max],...
set(gca,'XTickLabel',[]);




text(15, 1000,['$\tau$ (mean)    = ~' num2str(nanmean(tau(:,SCENARIO)),2),' years'],'FontSize',18,'Interpreter','latex');
text(15, 900,['$\tau$ (median) = ~' num2str(nanmedian(tau(:,SCENARIO)),2),' years'],'FontSize',18,'Interpreter','latex');
text(15, 800,['$\tau$ (s.d.)       = ~' num2str(nanstd(tau(:,SCENARIO)),2),' years'],'FontSize',18,'Interpreter','latex');

tau=tau_2x_2050;

title('Histogram plots of \tau','FontSize',18,'FontName',FontName);
%Subplot D
subplot('position',[xp1+width+xsep yp1 width2 height]); hold on; box on;

ybar=NaN(length(xe),4);

% cmap2=cmap([3 4 2 1],:);
cmap2=cmap([2 3 1 4],:);

for i=1:length(xe)-1
   ybar(i,1)=sum(idx==2 & tau(:,SCENARIO)>=xe(i) & tau(:,SCENARIO)<xe(i+1));
   ybar(i,2)=sum(idx==3 & tau(:,SCENARIO)>=xe(i) & tau(:,SCENARIO)<xe(i+1));
   ybar(i,3)=sum(idx==1 & tau(:,SCENARIO)>=xe(i) & tau(:,SCENARIO)<xe(i+1));
   ybar(i,4)=sum(idx==4 & tau(:,SCENARIO)>=xe(i) & tau(:,SCENARIO)<xe(i+1));
end

ybar(end,1)=sum(idx==2 & tau(:,SCENARIO)>xe(end));
ybar(end,2)=sum(idx==3 & tau(:,SCENARIO)>xe(end));
ybar(end,3)=sum(idx==1 & tau(:,SCENARIO)>xe(end));
ybar(end,4)=sum(idx==4 & tau(:,SCENARIO)>xe(end));

bh=bar(xe,ybar,'stacked'); 
% colormap(gca,cmap2);  
for i=1:4
    bh(i).FaceColor = cmap2(i,:);
end 

xlabel('\tau [years]','FontSize',14);
ylabel('Count (No. of stations)','FontSize',14);

set(gca,'FontSize',14,...
    'XLim',[5 35],'XTick',[0:5:50],...
    'YLim',[0 1250],...
    'FontName',FontName,...
    'Layer','top');
% 'YLim',[0 count_max],...

annotation(gcf,'textbox',...
    [0.501 0.923571428571428 0.0879999999999999 0.075],...
    'String',{'\tau [years]'},...
    'LineStyle','none',...
    'FontSize',18,...
    'FontName',FontName,...
    'FitBoxToText','off');

annotation(gcf,'textbox',...
    [0.0790000000000001 0.541428571428571 0.0359999999999999 0.0557142857142853],...
    'VerticalAlignment','middle',...
    'String',{'A'},...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',30,...
    'FontName','Nomada Incise',...
    'FitBoxToText','off',...
    'BackgroundColor',[0.8 0.8 0.8]);
annotation(gcf,'textbox',...
    [0.0790000000000001 0.0914285714285713 0.0359999999999999 0.0557142857142853],...
    'VerticalAlignment','middle',...
    'String','C',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',30,...
    'FontName','Nomada Incise',...
    'FitBoxToText','off',...
    'BackgroundColor',[0.8 0.8 0.8]);
annotation(gcf,'textbox',...
    [0.658000000000001+0.267 0.882857142857143 0.0359999999999999 0.0557142857142854],...
    'VerticalAlignment','middle',...
    'String','B',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',30,...
    'FontName','Nomada Incise',...
    'FitBoxToText','off',...
    'BackgroundColor',[0.8 0.8 0.8]);
annotation(gcf,'textbox',...
    [0.658000000000001+0.267  0.432857142857143 0.0359999999999999 0.0557142857142855],...
    'VerticalAlignment','middle',...
    'String','D',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',30,...
    'FontName','Nomada Incise',...
    'FitBoxToText','off',...
    'BackgroundColor',[0.8 0.8 0.8]);
text(15,1000,['$\tau$ (mean)    = ~' num2str(nanmean(tau(:,SCENARIO)),2),' years'],'FontSize',18,'Interpreter','latex');
text(15,900,['$\tau$ (median) = ~' num2str(nanmedian(tau(:,SCENARIO)),2),' years'],'FontSize',18,'Interpreter','latex');
text(15,800,['$\tau$ (s.d.)       = ~' num2str(nanstd(tau(:,SCENARIO)),2),' years'],'FontSize',18,'Interpreter','latex');

% print(gcf,'figures/Figure7_TIME_SCALE.jpg','-djpeg','-r300');
exportgraphics(gcf,'figures/Figure8.png','Resolution',300)


%check

% figure; hold on; box on;
% 
% id=lat>50;
% 
% geoshow('landareas.shp', 'FaceColor',[0.6 0.6 0.6],'EdgeColor',[0.6 0.6 0.6]); 
% axis equal; axis tight;
% scatter(lon(id),lat(id),30*ones(size(lat(id))),tau(id,SCENARIO),'filled'); colormap(gca,flipud(jet)); caxis([0 25]); colorbar;
% axis([-180 -60 0 80]);
% set(gca,'FontSize',14,'XTick',[-180:20:-60],'YTick',[0:10:80],'Layer','top');
% 
% tau(id,1);


%% Exceedance Plot with Time ODDS - (Figure 6)

% quick check of SLR curves

han1=figure(6666); set(han1,'PaperpositionMode','auto','Position',[200 50 800 1000],'color','w');

TEMP_PRJ=[Tmin, Tmean, Tmax]; % ????
% TIME_min=datenum(2000,1,1); % Global 
TIME_min=T_proj_local.time_start; % Local
TIME_max=datenum(2100,1,1);
TIME=(TIME_min:100:TIME_max)';

for SCENARIO=1:3 % 1=RCP 8.5, 2=RCP 4.5, 3=RCP 2.6
    
       
    for i=1:1:len
        
        
        id=ID(i); % get the specific station
        
%         [SLR_min(i),tid]=min(S(id).SLR_med(:,SCENARIO));
        [TEMP_min(i),tid]=min(TEMP_PRJ(:,SCENARIO));
        
        %plot(S(id).t_SLR,S(id).SLR_med(:,SCENARIO),'b-',S(id).t_SLR(tid),S(id).SLR_med(tid,SCENARIO),'ro'); datetick('x');
    end
    
    id_neg=find(TEMP_min<-0.01);
    
    % S(ID(id_neg)).name

    E_TIME=NaN(length(TIME),len);
    
    E_TIME_DATA=NaN(length(TIME),len);
    
    E_TIME_GEV=NaN(length(TIME),len);
    
    ODDS_TIME=NaN(length(TIME),len);
    ODDS_TIME_DATA=NaN(length(TIME),len);
    ODDS_TIME_GEV=NaN(length(TIME),len);
    
    TR=50;
    
    x1=(0:0.001:2*max(x_data))';
    
%     for i=1:len
    for i=1:1:len  % In this case not all the stations are being considered
        [SCENARIO i]
        id=ID(i); % get the specific station
        
        E0=S(id).ri/TR; % get the exceedance probability of the TR-year event
        
        % GEV
        mu=S(id).mu;
        sig=S(id).sigma;
        k=S(id).k;

%         SLR=interp1(S(id).t_SLR,S(id).SLR_med(:,SCENARIO),TIME);
%         GW=interp1(ti,TEMP_PRJ(:,SCENARIO),TIME); % Global
        GW=interp1(S(id).ti,S(id).TR_scenario(:,SCENARIO),TIME);
        
        % empirical exceedace probability of all data
        [Etmp,ia,~] = unique(S(id).Exceedance,'last');
        xtmp=S(id).x(ia);
        
        x0=x_Tr(TR,mu,sig,k); % find the water level corresponding to the 50-year event
        
        Ei=interp1(xtmp,Etmp,x0,'linear');
        
        for kk=1:length(GW)
            
            % find empirical exceedance probability increase of the all values
            E_TIME(kk,i)=interp1(xtmp+GW(kk),Etmp,x0);
                    
        end
        
        ODDS_TIME(:,i)=E_TIME(:,i)./(1-E_TIME(:,i));
        
        ODDSi=Ei./(1-Ei);
        
        color1=cmap(idx(i),:); % color based on clusters
        
        if ismember(i,id_neg)
            color1='k';           % single color
        end
        
        ALPHA=0.15;
        
        xp0=0.1;
        yp0=0.66;
        width=0.36;
        height=0.27;
        xsep1=0.125;
        xsep2=0.02;
        ysep=0.03;
        
        % This can be taken out from the for loop since I don't have a
        % separate GW for each station
        subplot('position',[xp0 yp0-(SCENARIO-1)*(height+ysep) width height]); hold on; box on;
        han1=plot(TIME,GW,'Color',color1); han1.Color(4)=ALPHA;  % set transparency
%         han1=plot(TIME,GW,'Color','k');
        set(han1,'LineWidth',2);


%         if SCENARIO==1;
%             title({'ensemble mean GW'});
%             ylabel('RCP 8.5 SLR [m]','FontSize',14);
%         elseif SCENARIO==2;
%             ylabel('RCP 4.5 SLR [m]','FontSize',14);
%         elseif SCENARIO==3;
%             ylabel('RCP 2.6 SLR [m]','FontSize',14);  
%         end
        if SCENARIO==1
            title({'Local temperature projections'},'FontName',FontName,'FontSize',18);
            ylabel('Low-end RCP 8.5 TR [\circC]','FontSize',14,'FontName',FontName);
        elseif SCENARIO==2
            ylabel('Intermediate RCP 8.5 TR [\circC]','FontSize',14,'FontName',FontName);
        elseif SCENARIO==3
            ylabel('High-end RCP 8.5 TR [\circC]','FontSize',14,'FontName',FontName);
        end
        
%         axis([TIME_min TIME_max -0.5 2.5]);
        axis([TIME_min TIME_max -0.5 10.5]);
        
        set(gca,'FontSize',14,'FontName',FontName);
        xlabel('Year','FontSize',14,'FontName',FontName);
        datetick('x');
        if SCENARIO<3 
            set(gca,'xticklabel',[]);
            xlabel('');
        end

        if i == i(end)
            if SCENARIO == 1
                AVG_Global = Tmin(1:1212,:);
                AVG_Local = mean(T_low_stations(:,1:101),1);
                
            elseif SCENARIO == 2
                AVG_Global = Tmean(1:1212,:);
                AVG_Local = mean(T_mean_stations(:,1:101),1);
            elseif SCENARIO == 3
                AVG_Global = Tmax(1:1212,:);
                AVG_Local = mean(T_high_stations(:,1:101),1);
            end

            AVG_Global_int = interp1(1:length(AVG_Global),AVG_Global',1:length(AVG_Global)/length(TIME):length(AVG_Global));
            AVG_Local_int = interp1(1:length(AVG_Local),AVG_Local,linspace(1,length(AVG_Local),length(TIME)));

            han11 = plot(TIME,AVG_Global_int,'--k','LineWidth',2);
            han22 = plot(TIME,AVG_Local_int,'--b','LineWidth',2);

            if SCENARIO == 1
                legend([han11 han22],'Global average','North America average','FontSize',13,'FontName',FontName);
            end

        end


        subplot('position',[xp0+xsep1+1*(width+xsep2) yp0-(SCENARIO-1)*(height+ysep) width height]); hold on; box on;
        han1=plot(TIME,ODDS_TIME(:,i)/ODDSi,'Color',color1); han1.Color(4)=ALPHA;  % set transparency
        set(han1,'LineWidth',2);
        if SCENARIO==1
            title({'All data'},'FontName',FontName,'FontSize',18);
        end
        set(gca,'FontSize',14,'FontName',FontName);
        datetick('x');
        xlabel('Year','FontSize',14,'FontName',FontName);
        if SCENARIO<3   
            set(gca,'xticklabel',[]);
            xlabel('');
        end
        
    end  


%     TIME_min=datenum(2000,1,1);
%     TIME_max=datenum(2100,1,1);
    
    subplot('position',[xp0+xsep1+width+xsep2 yp0-(SCENARIO-1)*(height+ysep) width height]); hold on; box on;

    text(731306.861571456, 82.6691307776895,'1 year','Rotation',88,'FontSize',14,'color','m','FontName',FontName);
    text(739435.75474142, 59.3537289020805 ,'5 years','Rotation',67,'FontSize',14,'color','b','FontName',FontName);
    text(746790.844100608, 33.7572473193161 ,'10 years','Rotation',50,'FontSize',14,'color','k','FontName',FontName);
    text(755415.747120308, 8.80139136611196 ,'25 years','Rotation',25,'FontSize',14,'color',[0.4940 0.1840 0.5560],'FontName',FontName);
    text(757698.809620308, 3.64269028542452 ,'50 years','Rotation',15,'FontSize',14,'color','r','FontName',FontName);
       
    plot(TIME,2.^((TIME-TIME_min)/(5*365.25)),'--','LineWidth',3,'Color','b'); hold on
    plot(TIME,2.^((TIME-TIME_min)/365.25),'--','LineWidth',3,'Color','m')
    plot(TIME,2.^((TIME-TIME_min)/(25*365.25)),'--','LineWidth',3,'Color',[0.4940 0.1840 0.5560])
    plot(TIME,2.^((TIME-TIME_min)/(50*365.25)),'--','LineWidth',3,'Color','r')
    plot(TIME,2.^((TIME-TIME_min)/(10*365.25)),'--','LineWidth',3,'Color','k')
    
    ylabel('$\frac{\mathrm{Odds}}{\mathrm{Odds}_0}$','FontSize',30,'Interpreter','latex');
    
    set(gca,'Yscale','log','ytick',[1 2 4 8 16 32 64 128 256],'YMinorTick','off');
    
    axis([TIME_min-185 TIME_max 1 256]);

end

annotation(gcf,'textbox',...
    [0.08 0.926469945355196 0.908333333333336 0.063333333333335],...
    'String',{'Normalized increase in the odds of exceedance (i.e., "doubling")'},...
    'LineStyle','none',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',18,...
    'FontName',FontName,...
    'FitBoxToText','off');

SCENARIO=1;
xp0=0.10875;  yp0=0.883926229508197;
xsep1=0.425;
ysep1=0.22;

annotation(gcf,'textbox',[xp0                       yp0-(SCENARIO-1)*(height+ysep) 0.04375 0.04],'String',{'A'},'FontWeight','bold','FontSize',30,'FitBoxToText','off','BackgroundColor',[0.8 0.8 0.8],'HorizontalAlignment','center','VerticalAlignment','middle','FontName','Nomada Incise');
annotation(gcf,'textbox',[xp0+xsep1+1*(width+xsep2) yp0-ysep1-(SCENARIO-1)*(height+ysep) 0.04375 0.04],'String',{'B'},'FontWeight','bold','FontSize',30,'FitBoxToText','off','BackgroundColor',[0.8 0.8 0.8],'HorizontalAlignment','center','VerticalAlignment','middle','FontName','Nomada Incise');

SCENARIO=2;
annotation(gcf,'textbox',[xp0                       yp0-(SCENARIO-1)*(height+ysep) 0.04375 0.04],'String',{'C'},'FontWeight','bold','FontSize',30,'FitBoxToText','off','BackgroundColor',[0.8 0.8 0.8],'HorizontalAlignment','center','VerticalAlignment','middle','FontName','Nomada Incise');
annotation(gcf,'textbox',[xp0+xsep1+1*(width+xsep2) yp0-ysep1-(SCENARIO-1)*(height+ysep) 0.04375 0.04],'String',{'D'},'FontWeight','bold','FontSize',30,'FitBoxToText','off','BackgroundColor',[0.8 0.8 0.8],'HorizontalAlignment','center','VerticalAlignment','middle','FontName','Nomada Incise');

SCENARIO=3;
annotation(gcf,'textbox',[xp0                       yp0-(SCENARIO-1)*(height+ysep) 0.04375 0.04],'String',{'E'},'FontWeight','bold','FontSize',30,'FitBoxToText','off','BackgroundColor',[0.8 0.8 0.8],'HorizontalAlignment','center','VerticalAlignment','middle','FontName','Nomada Incise');
annotation(gcf,'textbox',[xp0+xsep1+1*(width+xsep2) yp0-ysep1-(SCENARIO-1)*(height+ysep) 0.04375 0.04],'String',{'F'},'FontWeight','bold','FontSize',30,'FitBoxToText','off','BackgroundColor',[0.8 0.8 0.8],'HorizontalAlignment','center','VerticalAlignment','middle','FontName','Nomada Incise');

% print(gcf,'figures/Figure6_rate_of_increase_with_time_ODDS.jpg','-djpeg','-r300');
exportgraphics(gcf,'figures/Figure7.png','Resolution',300)

%% Appendix 5
id=~isnan(sigma_2x);
[S(ID(id)).sigma];
sigma_2x(id);

figure(7); set(gcf,'PaperpositionMode','auto','Position',[200 100 500 500],'color','w'); hold on; box on;
scatter([S(ID(id)).sigma], sigma_2x(id),55 ,'k','filled' );
scatter([S(ID(id)).sigma], sigma_2x(id),40 ,idx(id),'filled' );
colormap(cmap);
set(gca,'FontSize',16,'FontName',FontName);
axis equal;
axis([0,4,0,4])
yticks([0 1 2 3 4])

han1=plot([0 4],[0 4],'b--'); set(han1,'LineWidth',2);


xlabel(append('$\sigma$', ' [$^\circ$C]'), 'Interpreter', 'LaTeX', 'FontSize', 18);
ylabel(append('$\tilde{\sigma}$',' [$^\circ$C]'), 'Interpreter', 'LaTeX', 'FontSize', 18);

RR = corr([S(ID(id)).sigma]',sigma_2x(id)); RR = num2str(RR,'%.2f');

text(3.06242268041237, 0.149484536082474,append('$R$ = ',RR),'FontSize',20,'interpreter','latex');


% print(gcf,'figures/Appendix5.jpg','-djpeg','-r300');
exportgraphics(gcf,'figures/FigureA6.png','Resolution',300)

%% save all variables in a .mat file

save('All_Variables_nosmoothing.mat','-v7.3','-nocompression');

%% Extra figures for reviewers (not used in manuscript)



