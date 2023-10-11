function S=conv_stations2S(stations,states)
fields=fieldnames(stations(1).list);
c = cell(length(fields),sum([stations.numberOfStations]));
S = cell2struct(c,fields);
for i = 1:numel(stations)
    for ii = 1:numel(stations(i).list)
        stationnum=sum([states(1:i-1).numberOfStations])+ii;
        try
            S(stationnum)=stations(i).list(ii);
        catch
            fields=fieldnames(stations(i).list(ii));
            for j=1:numel(fields)
                S(stationnum).(fields{j})=stations(i).list(ii).(fields{j});
            end
        end
    end
end