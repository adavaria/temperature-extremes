function [stations, states, countries]=checking_stations_offline(country,element)
% data_path='/Volumes/ALI - WD/Dr. Vitousek/data/ghcn/';
data_path = 'data/ghcnd/';
if ispc
    data_path='C:\Users\adava\Downloads\ghcnd_all\';
end
if nargin < 1
    element='TMAX';
    %country={'AE','AG','AF'};
    country='US';
end
if ~iscell(country)
    country={country};
end
if size(country,2)>1
    country=country';
end

filename='ghcnd-countries.txt';
fileID = fopen([data_path filename],'r');
st=strtrim(textscan(fileID,'%2s%s','whitespace','\n'));
fclose(fileID);
countries=table(st{1},st{2});
countries.Properties.VariableNames={'country_code','country_name'};
if strcmpi(country,'all')
    country=countries.country_code;
end


filename='ghcnd-stations.txt';
fileID = fopen([data_path filename],'r');
if fileID<0
    error('Couldn''t locate the file');
end
t=textscan(fileID,'%11s%1s%8f%1s%9f%1s%6f%1s%2s%1s%30s%1s%3s%1s%3s%1s%5s','whitespace','');
fclose(fileID);
t{2}=extractBefore(t{1},3); % extracting country codes
s=cell2mat([t{1}]);
loc_t=[];
for i=1:numel(country)
    loc_t=cat(1,loc_t,find(strcmp(cellstr(s(:,1:2)),country{i}))); %within the country(ies)
end
T_t=table(t{1}(loc_t),strtrim(t{11}(loc_t)),num2cell(t{3}(loc_t)),...
    num2cell(t{5}(loc_t)),num2cell(t{7}(loc_t)),t{2}(loc_t),t{9}(loc_t));
% string() can be used if I didn't like to have them as cell in my table
T_t.Properties.VariableNames={'id','name','lat','lon','ele','country_code','st'};
nonst='N/A';
T_t.st(strcmp(T_t.st,'  '))={nonst};


filename='ghcnd-inventory.txt';
fileID = fopen([data_path filename],'r');
in=textscan(fileID,'%11s%20s%4s%1s%4d%1s%4d','whitespace','');
fclose(fileID);
s=cell2mat([in{1}]);
loc_cont=[];
for i=1:numel(country)
    loc_cont=cat(1,loc_cont,find(strcmp(cellstr(s(:,1:2)),country{i}))); %within the country
end
loc_el=find(strcmp(in{3},element)); % Element restriction
loc_per=find(in{7}-in{5}>=10); % Period restriction (years)
loc_in=intersect(intersect(loc_el,loc_per),loc_cont); % mixing constraints
T_in=table(in{1}(loc_in),double(in{5}(loc_in)),double(in{7}(loc_in))...
    ,double(in{7}(loc_in)-in{5}(loc_in)),cell(size(loc_in)));
T_in.Properties.VariableNames={'id','firstyear','lastyear','period','S'};

filename='ghcnd-states.txt';
fileID = fopen([data_path filename],'r');
st=strtrim(textscan(fileID,'%2s%s','whitespace','\n'));
fclose(fileID);
st{1}{end+1}=nonst;st{2}{end+1}=nonst;
states=table(st{1},st{2});
states.Properties.VariableNames={'st','state'};
states=sortrows(states,'state');

stations=innerjoin(T_t,T_in);
stations=innerjoin(stations,states);
tmp=categorical(stations.st);
states=innerjoin(table(cellstr(unique(tmp)),'VariableNames',{'st'}),states);
states=[states table(countcats(tmp),'VariableNames',{'numberOfStations'})];


countries=innerjoin(countries,table(country,'VariableNames',{'country_code'}));
tmp=categorical(stations.country_code);
countries=innerjoin(countries,...
    table(unique(stations.country_code),'VariableNames',{'country_code'})); % To eliminate countries with no stations
countries=[countries table(countcats(tmp),'VariableNames',{'numberOfStations'})];

stations=innerjoin(stations,countries);

% save stations.mat stations states

% all=struct('id',t{1}(loc_us),'name',t{11}(loc_us),'lat',num2cell(t{3}(loc_us)),'lan',num2cell(t{5}(loc_us)),...
%     'ele',num2cell(t{7}(loc_us)),'state',t{9}(loc_us));
% disp(all);