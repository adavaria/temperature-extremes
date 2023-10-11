clearvars -except saved_stations states
clc
if ~exist('saved_stations','var')||~exist('states','var')
    load selected_stations_completed.mat;
end
current_path=pwd;
irun=5;
data_path='/Volumes/ALI - WD/Dr. Vitousek/data/ghcn/ghcnd_all/ghcnd_all/';
str='/';
if ispc
    str='\';
    data_path='C:\Users\adava\Downloads\ghcnd_all\ghcnd_all\';
    
end
%cd(data_path)
%fileList=dir;
%cd ..
save_path=[data_path(1:end-10) 'run' num2str(irun) str];
% winopen(save_path)
%if exist(['run' num2str(irun)],'dir')==0
if exist(save_path,'dir')==0
    mkdir([save_path 'run' num2str(irun)]);
else
%     str=input('Do you want to delete previous results?, Y/N [Y]:','s');
%     if upper(str)=='Y'
%         delete([save_path '*']);
%     end
end
%cd(current_path);



total_time=0;
stationnum=0;

cont=0;
if cont==1
    load([data_path 'get_data_tmp.mat'])
else
    i=1;
    ii=1;
end


time_total0=cputime;
%while i<=3 % i refers to states
while i<=numel(states) % i refers to states
    fprintf('%d-%s (%s): Geting data...\n',i,states(i).name,states(i).stateABR);
    tic
    time_state0=cputime;
    %state=stations(i);
    %for ii=1:3
    for ii=1:numel(stations(i).list)
        stationnum=sum([states(1:i-1).numberOfStations])+ii;
        station=stations(i).list(ii);
        STATION=station.id;
        fprintf('\t%d- %s, %s: ',ii,STATION,station.name)
        %fileID=find(strcmp({fileList.name},[STATION(7:17) '.dly']), 1);
        %if isempty(fileID)
        if exist([data_path [STATION(7:17) '.dly']],'file')<1
            warning('Couldnt find the datafile');
            stations(i).list(ii).QC='Couldn''t find the datafile';
            stations(i).list(ii).completed=false; %#ok<*SAGROW>
            fprintf(2,'Error: %s\n',stations(i).list(ii).QC)
            continue
        end
        [h_all, t_all, flags]=read_dly([STATION(7:17) '.dly'],'TMAX',data_path);
        % Checing the coverage and period
        if (numel(h_all)-sum(sum(isnan(h_all))))/(t_all(end)-t_all(1))<.8
            stations(i).list(ii).QC='No Enough Coverage';
            stations(i).list(ii).completed=false; %#ok<*SAGROW>
            fprintf(2,'Error: %s\n',stations(i).list(ii).QC)
            continue
        end
        period=str2num(datestr(t_all(end),'yyyy'))-str2num(datestr(t_all(1),'yyyy')); %#ok<ST2NM>
        if period<10
            stations(i).list(ii).QC='Less than 10 years data';
            stations(i).list(ii).completed=false; %#ok<*SAGROW>
            fprintf(2,'Error: %s (%d)\n',stations(i).list(ii).QC,period)
            continue
        end
        
        h_all=h_all/10; % Unit Converstion from 10th of C to C
        % finding outliers and elemination
        %tmp=find(abs(T_all)>100);
        tmp=find(~isspace(flags.qflags));
        h_outliers=h_all(tmp);
        t_outliers=t_all(tmp);
        
        % Normal Outlier
%         OL=isoutlier(T_all);
%         OL_loc=find(OL);
%         plot(t_all(OL_loc),T_all(OL_loc),'r*');
%         hold on
%         Tp=T_all';
%         tp=t_all';
%         f=find(diff(tp)>1);
%         tp(f)=nan;Tp(f)=nan;
%         plot(tp,Tp,'-b.');
%         figure
%         plot(t_all(tmp),T_all(tmp),'g*')
%         hold on
%         plot(t_all,T_all,'-b.');
        
        
        h_all(tmp)=NaN;%t_all(tmp)=NaN;
        % selecting top events
        [h_events, t_events]=top_events(h_all,t_all,20,3);
        
        % Processing
        Nbins=500;
        [hist_N,hist_edges] = histcounts(h_all,Nbins);
        [F,edges] = histcounts(h_all,Nbins,'Normalization','cdf');
        [Ystart,~,~]=datevec(datenum(station.mindate));
        [Ystop,~,~]=datevec(datenum(station.maxdate));
        years=Ystart:Ystop;
        years_no_data=setdiff(Ystart:Ystop,str2num(datestr(t_events(:,1),'yyyy'))); %#ok<ST2NM>
        avg=nanmean(h_all(:));
        summer_avg=nan(numel(years),1);
        winter_avg=nan(numel(years),1);
        [yy,mm,~]=datevec(t_all);
        for j=1:numel(years)
            ind=(mm==12|mm==1|mm==2)&yy==years(j);
            if ind>2*30
                winter_avg(j)=nanmean(h_all(ind));
            end
            ind=and(or(mm==6,or(mm==7,mm==8)),yy==years(j));
            if ind>2*30
                summer_avg(j)=nanmean(h_all(ind));
            end
        end
        h_coverage=1-sum(sum(isnan(h_all)))/numel(h_all);
        str_fig=sprintf('%s%d-%d_%s.png',save_path,i,ii,station.name);
        T_gev=h_events(:,1:3);
        [parmhat, parmci] = gevfit(T_gev(:));
%         k=parmhat(1);sigma=parmhat(2);mu=parmhat(3);
%         k_ci=parmci(:,1);sigma_ci=parmci(:,2);mu_ci=parmci(:,3);
        
%         k=nan(size(T,1),1);sigma=k;mu=k;
%         warning ('off','all');
%         for j=1:size(T,1)
%             parmhat = gevfit(T(j,1:5));
%             k(j)=parmhat(1);sigma(j)=parmhat(2);mu(j)=parmhat(3);
%         end
%         warning ('on','all');
        

        stations(i).list(ii).h_coverage=h_coverage;
        stations(i).list(ii).completed=true;
        stations(i).list(ii).QC='Everything seems ok.';
        stations(i).list(ii).t_events=t_events;
        stations(i).list(ii).h_events=h_events;
        stations(i).list(ii).years=years;
        stations(i).list(ii).Nyears=numel(years);
        stations(i).list(ii).Npeaks=size(h_events,2);
        stations(i).list(ii).yearsWithData=size(h_events,1);
        stations(i).list(ii).Npeaks=numel(t_events);
        stations(i).list(ii).years_no_data=years_no_data;
        stations(i).list(ii).edges=edges';
        stations(i).list(ii).x=0.5*(edges(2:end)+edges(1:end-1))';
        stations(i).list(ii).Exceedance=1-F';
        stations(i).list(ii).hist_N=hist_N';
        stations(i).list(ii).hist_edges=hist_edges';
        stations(i).list(ii).winter_average=winter_avg;
        stations(i).list(ii).summer_average=summer_avg;
        stations(i).list(ii).k=parmhat(1);
        stations(i).list(ii).sigma=parmhat(2);
        stations(i).list(ii).mu=parmhat(3);
        stations(i).list(ii).k_ci=parmci(:,1);
        stations(i).list(ii).sigma_ci=parmci(:,2);
        stations(i).list(ii).mu_ci=parmci(:,3);
        stations(i).list(ii).FIG=str_fig;
        
        plotting=1;
        if ~plotting
            t_events=toc;
            fprintf('Successful. (%d sec.)\n',floor(t_events));
            continue
        end
        
        fig=figure;
        fig.Position(1:4)=[200 150 800 800];
        subplot(3,1,1)
        t_all_min=min(t_all(:));
        t_all_max=max(t_all(:));
        selected_outliers=find(h_outliers>min(h_all(:)) ...
            & h_outliers<max(h_all(:)));
        T_nan=zeros(t_all_max-t_all_min+1,1)+avg;
        t_nan=t_all(1):t_all(end);
        [~,ia,~] = intersect(t_nan,t_all(~isnan(h_all)));
        T_nan(ia)=nan;
        Tp=h_all';
        tp=t_all';
        f=find(diff(tp)>1);
        tp(f)=nan;Tp(f)=nan;
        hold on
        plot(tp,Tp,'-b.',t_events,h_events,'go');
        plot(t_outliers(selected_outliers),h_outliers(selected_outliers),'r*');
        a2=axis;
        plot(t_all_min:t_all_max,T_nan,'r','LineWidth',10)
        datetick('x','yyyy','keeplimits');
        xlabel('Time (year)');ylabel('Max. Temp. (C)');
        %legend('All Data','T Events','Out-Liers','Missing Data')
        hold off
        
        subplot(3,1,2)
        Tp=h_events';
        tp=t_events';
        [~,ind]=sort(tp);
        ind=repmat(0:size(t_events,1)-1,[size(t_events,2),1])*size(t_events,2)+ind;
        plot(tp(ind),Tp(ind),'-*','LineWidth',1.5)
        
        % [~,ind]=sort(t(:));
        % plot(t(ind),T(ind),'LineWidth',1.5);
        hold on
        plot([datenum([years_no_data' ones(numel(years_no_data),2)])'...
            ;datenum([years_no_data' ones(numel(years_no_data),2)])'+365],...
            ones(2,numel(years_no_data))+mean(h_events(:)),'r','LineWidth',20);
        str=sprintf('Station ID: %s \n Station Name: %d-%d- %s',station.id,i,ii,station.name);
        title(str)
        hold off
        a1=axis;
        datetick('x','yyyy','keeplimits');
        xlabel('Time (year)');ylabel('Max. Temp. (C)');
        axis([a2(1) a2(2) a1(3) a1(4)])
        
        subplot(3,1,3);
        bar((hist_edges(1:end-1)+hist_edges(2:end))/2,hist_N,3)
        xlabel('Max. Temp. (C)')
        ylabel('Total Number of Events')
        
%         subplot(2,2,4)
%         str=evalc(sprintf('disp(stations(%d).list(%d))',i,ii));
%         f=find(double(str)==10);
%         text(0,.5,str(1:f(15)));
%         text(.6,0.5,str(f(15)+1:end))
%         axis off
        
        
        saveas(fig,str_fig);
        close(fig);
        ti=toc;
        fprintf('Successful. (%d sec.)\n',floor(ti));
    end
    time_state0=floor(cputime-time_state0);
    fprintf('%s (%s): Completed (%d min., %d sec.)\n',states(i).name,...
        states(i).stateABR,floor(time_state0/60),mod(time_state0,60));
    state=states(i);
    state.list=stations(i).list;
    save([save_path 'get_data_tmp.mat'])
    save([save_path state.name '.mat'],'state')
    save([save_path 'stations.mat'],'stations')
    fprintf('\n')
    %ii=1;
    i=i+1;
end
disp('Converting stations into S structure');
S=conv_stations2S(stations,states);
save([save_path 'S.mat'],'S')
save('S.mat','S');

time_total0=cputime-time_total0;
fprintf('Process completed. Total time = %.0f hr., %.0f min.\n',time_total0/3600,...
    (-floor(time_total0/3600)*3600+time_total0)/60);
disp(['Click ',...
        '<a href = "matlab: [s,r] = system(''explorer ',save_path,' &'');">',...
        'here','</a>',...
        ' to open the results location.'])
    