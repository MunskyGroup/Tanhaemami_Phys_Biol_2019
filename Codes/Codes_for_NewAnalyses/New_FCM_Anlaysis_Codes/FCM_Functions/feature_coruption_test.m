for Target_Ind = 1:14; %Choose the target parmeter index(3 and 10 are targets).
    figure(Target_Ind); clf; %Create figure, clear it.
    for k=1:12 %For time points 1 to 12.
        try
            subplot(3,4,k)
            a=log(DATA.Stained(k).DATA(:,Target_Ind));
            ksdensity(a(a~=Inf & (a~=-Inf))); %Plot a kernel density of stained data.
            hold on
            b=log(DATA.Unstained(k).DATA(:,Target_Ind));
            ksdensity(b(b~=Inf & (b~=-Inf))); %Plot a kernel density of unstained data.
            legend('stained','unstained');
        catch
        end
    end
end