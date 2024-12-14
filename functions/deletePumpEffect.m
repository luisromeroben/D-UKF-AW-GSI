function out = deletePumpEffect(in,flowsRaw)
pressureRaw = in;
for s = 1:size(pressureRaw,2)
    if s >= 4
        pres = pressureRaw(:,s);
        f = flowsRaw(:,3);
        change = [];
        for i=2:length(pres)
            if f(i)-f(i-1)>30
                change = [change; i-1 i abs(pres(i)-pres(i-1)) 1];
            elseif f(i)-f(i-1)<-30
                change = [change; i-1 i abs(pres(i)-pres(i-1)) -1];
            end
        end
        
        pres(1:change(1,1)) = pres(1:change(1,1)) + change(1,3);
        
        for i=2:2:length(change)
            if i == length(change)
                pres(change(i,2):end) = pres(change(i,2):end) + change(end,3);
            else
                pres(change(i,2):change(i+1,1)) = pres(change(i,2):change(i+1,1)) + mean([change(i,3) change(i+1,3)]);
            end
        end
    else
        pres = pressureRaw(:,s);
    end
    
    out(:,s) = pres;
    
end

end