function fig=plot_station(data,STATION,element,data_path)
%% Prepating
% if ~exist('data','var')
%     data={};
% end
if isempty(data)
    if nargin<2
        STATION='USC00500172';
    end
    if nargin<3
        element='TMAX';
    end
    if nargin<4
        data_path='/Volumes/ALI - WD/Dr. Vitousek/data/ghcn/ghcnd_all/ghcnd_all/';
    %     str='/';
        if ispc
    %         str='\';
            data_path='C:\Users\adava\Downloads\ghcnd_all\ghcnd_all\';
        end
    end
    [h_all, t_all, flags]=read_dly([STATION '.dly'],element,data_path);
    h_all=h_all/10; % Unit Converstion from 10th of C to C
    % finding outliers and elemination
    tmp=find(~isspace(flags.qflags));
    h_outliers=h_all(tmp);
    t_outliers=t_all(tmp);
    h_all(tmp)=NaN;%t_all(tmp)=NaN;
    % selecting top events
    [h_events, t_events]=top_events(h_all,t_all,20,3);
else
    h_all=data{1};
    t_all=data{2};
    h_events=data{3};
    t_events=data{4};
    h_outliers=data{5};
    t_outliers=data{6};
end


%% Processing
Nbins=500;
[hist_N,hist_edges] = histcounts(h_all,Nbins);
tmp=datevec(t_all([1 numel(t_all)]));
Ystart=tmp(1);
Ystop=tmp(2);
years_no_data=setdiff(Ystart:Ystop,str2num(datestr(t_events(:,1),'yyyy'))); %#ok<ST2NM>
avg=nanmean(h_all(:));

%% Plotting
fig=figure;
fig.Position(1:4)=[200 150 800 800];
subplot(3,2,[1 2])
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
xlabel('Time (year)');ylabel('All Temp. (C)');
str=sprintf('Station ID: %s',STATION);
title(str)
hold off

subplot(3,2,[3 4])
Tp=h_events';
tp=t_events';
[~,ind]=sort(tp);
ind=repmat(0:size(t_events,1)-1,[size(t_events,2),1])*size(t_events,2)+ind;
plot(tp(ind),Tp(ind),'-*','LineWidth',1.5)
hold on
plot([datenum([years_no_data' ones(numel(years_no_data),2)])'...
    ;datenum([years_no_data' ones(numel(years_no_data),2)])'+365],...
    ones(2,numel(years_no_data))+mean(h_events(:)),'r','LineWidth',20);
hold off
a1=axis;
datetick('x','yyyy','keeplimits');
xlabel('Time (year)');ylabel('Max. Temp. (C)');
axis([a2(1) a2(2) a1(3) a1(4)])

subplot(3,2,5);
bar((hist_edges(1:end-1)+hist_edges(2:end))/2,hist_N,4)
xlabel('Temp. (C)')
ylabel('Total Number of Events')
subplot(3,2,6);
bar((hist_edges(1:end-1)+hist_edges(2:end))/2,cumsum(hist_N),1)
xlabel('Temp. (C)')
ylabel('Cummulative Total Number of Events')

%         subplot(2,2,4)
%         str=evalc(sprintf('disp(stations(%d).list(%d))',i,ii));
%         f=find(double(str)==10);
%         text(0,.5,str(1:f(15)));
%         text(.6,0.5,str(f(15)+1:end))
%         axis off
