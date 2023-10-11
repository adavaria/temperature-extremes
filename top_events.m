function [h_events, t_events]=top_events(H,t,N,m)
% H: events
% t: time series
% N: number of top events
% m: min. distance of the events

if nargin<3
    N=20;
end
if nargin<4
    m=3;
end

H=reshape(H',[],1);
t=reshape(t',[],1);
t(isnan(H))=[];
H(isnan(H))=[];

years=unique(str2num(datestr(t,'yyyy'))); %#ok<ST2NM>
t_events=NaN(numel(years),N);
h_events=NaN(numel(years),N);

for Y=1:numel(years)
    START=datenum([num2str(years(Y)) '-01-01']);
    END=datenum([num2str(years(Y)) '-12-31']);
    ind=find(t>=START & t<=END);
    
    START_SUMMER=datenum([num2str(years(Y)) '-06-01']);
    END_SUMMER=datenum([num2str(years(Y)) '-08-31']);
    ind_summer=find(t>=START_SUMMER & t<=END_SUMMER);
    
    % Restrictions: period of 9 month and coverage of 6 month at
    % least and %80 of summer coverage
    
    if t(ind(end))-t(ind(1))<9*31 ||...
            numel(ind)-sum(isnan(H(ind)))<6*31 ||...
            numel(ind_summer)-sum(isnan(H(ind_summer)))<3*31*0.8
        continue
    end
    
    t_year=t(ind);
    T_year=H(ind);

    [~, locs]=findpeaks(T_year,'MinPeakDistance',m,'Npeaks',N,'SortStr','descend');   
    t_events(Y,1:length(t_year(locs)))=t_year(locs);
    h_events(Y,1:length(t_year(locs)))=T_year(locs);
end

t_events(sum(isnan(h_events),2)==N,:)=[];
h_events(sum(isnan(h_events),2)==N,:)=[];

%str2num(datestr(t,'yyyy')) %#ok<ST2NM>
% [Y,M,D]=datevec(t);
% disp([Y M D T*10])

% figure
% plot(t,T,'-b.',t_events(:,1:5),T_events(:,1:5),'ro')
% datetick('x','yyyy','keeplimits');

%years(sum(isnan(T_events),2)==N)=[];

% 
% figure
% plot(t(:),T(:));
% datetick('x','yyyy','keeplimits');