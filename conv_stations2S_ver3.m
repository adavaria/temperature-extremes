function S=conv_stations2S_ver3(stations)
S=repmat(struct,[size(stations,1) 1]);
for i=1:numel(S)
    try 
        S(i)=stations.S{i};
    catch
        tmp=stations.S{i};
        fields=fieldnames(tmp);
        for j=1:numel(fields)
            S(i).(fields{j})=tmp.(fields{j});
        end
    end
end