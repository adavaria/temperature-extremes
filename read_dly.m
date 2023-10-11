function [H,t,flags]=read_dly(filename,var,pathname)
if nargin==0
    filename='USC00025512.dly';
end
if nargin <2
    var='TMAX';
end
if nargin<3
    pathname='';
end
fileID = fopen([pathname filename],'r');
if fileID<1
    error('Couldn''t located the file.')
end

t=textscan(fileID,['%11s%4d%2d%4s' repmat('%5d%3s',[1,31])],'whitespace','');
dates=double([t{[2 3]}]);
n=size(dates,1);
dates=datenum([dates zeros(n,1)]);
dates=dates+repmat(1:31,[n 1]);
values=double([t{5:2:65}]);
values(values==-9999)=NaN;
flags_tmp=cell2mat([t{6:2:66}]);
var_pos=find(strcmp(t{4},var));
if isempty(var_pos)
    warning('No results found');
    H=[];
    t=[];
    return;
end
flags=struct;
flags.qflags=flags_tmp(var_pos,2:3:end);
H=values(var_pos,:);
t=dates(var_pos,:);

fclose(fileID);

% test dates:
% reshape(cellstr(datestr(t,'yyyy-mm-dd')),[size(t,1),31])
% plot(reshape(t',[],1),reshape(T',[],1),'-b.');
% plot(reshape(t_all',[],1),reshape(T_all',[],1),'-b.');
% test dates:
% reshape(cellstr(datestr(dates,'yyyy-mm-dd')),[n,31])