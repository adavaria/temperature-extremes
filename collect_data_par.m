%clearvars -except saved_stations states
clear
clc
close all

element='TMAX';
country={'US','MX','CA'};
i=1;
station0=1;
plotting=1;
cont=0;
irun=12;

% data_path='/Volumes/ALI - WD/Dr. Vitousek/data/ghcn/ghcnd_all/ghcnd_all/';
data_path = 'data/ghcnd/ghcnd_all/ghcnd_all/';
str='/';
if ispc
    str='\';
    %filesep
    data_path='C:\Users\adava\Downloads\ghcnd_all\ghcnd_all\';
end
save_path=[data_path(1:end-10) 'run_' num2str(irun) str];
% winopen(save_path)
if exist(save_path,'dir')==0
    mkdir(save_path);
%     else
%         str=input('Do you want to delete previous results?, Y/N [Y]:','s');
%         if upper(str)=='Y'
%             diary off
%             fclose('all');
%             delete([save_path '*']);
%         end
end

if cont==1
    load([save_path 'get_data_tmp.mat'])
else
    [stations, states, countries]=checking_stations_offline(country,element);
end

str='diary.txt';
if exist([save_path str],'file')~=0
    delete([save_path str]);
end
diary off
diary([save_path str]);
% p=gcp;
% disp(p);

time_total0=datetime;
BadChar = '<>:"/\|?*';
fprintf('Total number of stations: %d\n\n',size(stations,1));
fprintf('Countries: \n');
disp(countries)
N=100; % Batch size
Nmax=3; ri=1/Nmax;
total_batches=floor(size(stations,1)/N)+1;
while i<=total_batches
    time_batch0=datetime;
    list_ind=(i-1)*N+1:min(i*N,size(stations,1));
    list=table2struct(stations(list_ind,:));
    C=cell(numel(list_ind),1);
    parfor ii=station0:numel(list_ind)
%     for ii=station0:numel(list_ind)
        Ctmp=list(ii);
        STATION=list(ii).id;
        fprintf('\t%d (%d-%d)- %s, %s: ',ii+(i-1)*N,i,ii,STATION,list(ii).name)
        if exist([data_path [STATION '.dly']],'file')<1
            warning('Couldnt find the datafile');
            Ctmp.QC='Couldn''t find the datafile';
        end
        [h_all, t_all, flags]=read_dly([STATION '.dly'],element,data_path);
        % Checing the coverage and period
        if (numel(h_all)-sum(sum(isnan(h_all))))/(t_all(end)-t_all(1))<.8
            Ctmp.QC='No Enough Coverage';
        end
        period=str2num(datestr(t_all(end),'yyyy'))-str2num(datestr(t_all(1),'yyyy')); %#ok<ST2NM>
        if period<10
            Ctmp.QC='Less than 10 years data';
        end
        if isfield(Ctmp,'QC')
            Ctmp.completed=false; %#ok<*SAGROW>
            C{ii}=Ctmp;
            fprintf(2,'Error: %s\n',Ctmp.QC);
            continue;
        end
        h_all=h_all/10; % Unit Converstion from 10th of C to C
        % finding outliers and elemination
        tmp=find(~isspace(flags.qflags));
        h_outliers=h_all(tmp);
        t_outliers=t_all(tmp);   
        h_all(tmp)=NaN;%t_all(tmp)=NaN;
        % selecting top events
        [h_events, t_events]=top_events(h_all,t_all,20,3);
        if isempty(h_events)
            Ctmp.QC='No top-events constraints wes met';
            Ctmp.completed=false; %#ok<*SAGROW>
            C{ii}=Ctmp;
            fprintf(2,'Error: %s\n',Ctmp.QC);
            continue;
        end
        
        % Processing
        Nbins=500;
        [hist_N,hist_edges] = histcounts(h_all(~isnan(h_all)),Nbins);
        [F,edges] = histcounts(h_all(~isnan(h_all)),Nbins,'Normalization','cdf');
        %[Ystart,~,~]=datevec(datenum(list(ii).mindate));
        %[Ystop,~,~]=datevec(datenum(list(ii).maxdate));
        Ystart=list(ii).firstyear;
        Ystop=list(ii).lastyear;
        years=Ystart:Ystop;
        years_no_data=setdiff(Ystart:Ystop,str2num(datestr(t_events(:,1),'yyyy'))); %#ok<ST2NM>
        avg=nanmean(h_all(:));
        summer_avg=cat(2,nan(numel(years),1),years');
        winter_avg=cat(2,nan(numel(years),1),years');
        [yy,mm,~]=datevec(t_all);
        for j=1:numel(years)
            ind=(mm==12|mm==1|mm==2)&yy==years(j);
            if sum(ind(:))>3*30*0.8
                winter_avg(j,1)=nanmean(h_all(ind));
            end
            ind=and(or(mm==6,or(mm==7,mm==8)),yy==years(j));
            if sum(ind(:))>3*30*0.8
                summer_avg(j,1)=nanmean(h_all(ind));
            end
        end
        summer_avg=round(summer_avg,2);
        winter_avg=round(winter_avg,2);
        h_coverage=1-sum(sum(isnan(h_all)))/numel(h_all);
        str_fig=sprintf('%s%d-%d_%s.png',save_path,i,ii,list(ii).name);
        tmp=str_fig(numel(save_path)+1:end); % fixing bad charachters for saving
        tmp(ismember(str_fig(numel(save_path)+1:end),BadChar))='_';
        str_fig(numel(save_path)+1:end)=tmp;
        
        t_max=t_events(:,1:Nmax);t_max=t_max(:);
        h_max=h_events(:,1:Nmax);h_max=h_max(:);
        [h_max,id]=sort(h_max,'descend'); 
        t_max=t_max(id);
        
        T_gev=h_events(:,1:ri^(-1));
        lastwarn('');
        [parmhat, parmci] = gevfit(T_gev(:));
        if ~isempty(lastwarn)
            Ctmp.QC='GEV fit failed';
        end
        
        
        Ctmp.h_coverage=h_coverage;
        if isfield(Ctmp,'QC')
            Ctmp.completed=false;
        else
            Ctmp.completed=true;
            Ctmp.QC='Everything seems ok.';
        end
        
        Ctmp.t_events=t_events;
        Ctmp.h_events=h_events;
        Ctmp.years=years;
        Ctmp.Nyears=numel(years);
        Ctmp.Npeaks=size(h_events,2);
        Ctmp.yearsWithData=size(h_events,1);
        Ctmp.Npeaks=numel(t_events);
        Ctmp.years_no_data=years_no_data;
        Ctmp.edges=edges';
        Ctmp.x=0.5*(edges(2:end)+edges(1:end-1))';
        Ctmp.Exceedance=1-F';
        Ctmp.hist_N=hist_N';
        Ctmp.hist_edges=hist_edges';
        Ctmp.winter_avg_vec=winter_avg;
        Ctmp.summer_avg_vec=summer_avg;
        Ctmp.winter_avg=round(nanmean(winter_avg(:,1)),2);
        Ctmp.summer_avg=round(nanmean(summer_avg(:,1)),2);
        Ctmp.k=parmhat(1);
        Ctmp.sigma=parmhat(2);
        Ctmp.mu=parmhat(3);
        Ctmp.k_ci=parmci(:,1);
        Ctmp.sigma_ci=parmci(:,2);
        Ctmp.mu_ci=parmci(:,3);
        Ctmp.ri=ri;
        Ctmp.Nmax=Nmax;
        Ctmp.h_max=h_max;
        Ctmp.t_max=t_max;
        
        
        
        if plotting==1
            Ctmp.FIG=str_fig;
        else
            Ctmp.FIG='No plotting';
        end
        C{ii}=Ctmp;
        
        if plotting==1
            fig=plot_station({h_all,t_all,h_events,t_events,h_outliers,t_outliers},STATION);
            saveas(fig,str_fig);
            close(fig);
        end
        
        fprintf('Successful.\n');
    end
    stations(list_ind,:).S=C;
    time_batch0=floor(datetime-time_batch0);
    fprintf('Batch %d/%d Completed (%s.)\n',i,total_batches,time_batch0);
    fprintf('Estimated remaining time = %s\n',...
        (datetime-time_total0)/(i*N)*(max(size(stations,1)-i*N,0)));
    save([save_path 'get_data_tmp.mat'])
    fprintf('\n')
    i=i+1;
end
disp('Converting stations into S structure');
S=conv_stations2S_ver3(stations);
save([save_path 'stations.mat'],'stations')
save([save_path 'S' '.mat'],'S')
% save(['S_' country '.mat'],'S');

time_total0=datetime-time_total0;
fprintf('Process completed. Total time = %s \n',time_total0);
disp(['Click ',...
    '<a href = "matlab: [s,r] = system(''explorer ',save_path,' &'');">',...
    'here','</a>',...
    ' to open the results location.'])
diary off
delete(gcp('nocreate'));