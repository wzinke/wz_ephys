function ISIviewer(clustnum,handles)

global clustattrib;
global clustdata;

plotfig = figure('Color',get(0,'DefaultUicontrolBackgroundColor'), ... 
        'Tag','viewISIfig','NumberTitle','off','Name',['ISI viewer: Cluster ',num2str(clustnum)]); 
    
%h1 = axes('XScale','log');
h1 = axes('XScale','linear');
hold on
index = clustattrib.clusters{clustnum}.index;
times = clustdata.params(index,1);
times = times/(clustdata.UnitsPerSec);
ISI = diff(times);
ISI = ISI(find(ISI<.5));
ISI = [ISI;-ISI];
edges = [-.3:.0005:.3];
%dges = [[-.025:.0005:.025]];
N = histc(ISI,edges);


if (sum(N))
    axis([-.3 .3 0 max(N)]);
    phandle = bar(edges,N,'histc');
    set(phandle,'LineStyle','none');
    set(phandle,'FaceColor',[0 0 0]);
    xlabel('Time between events (seconds)');
    ylabel('Number of events');
    %set(phandle,'XScale','log');
else
    axis([-.3 .3 0 1]);
end