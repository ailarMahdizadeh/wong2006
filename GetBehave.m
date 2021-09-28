function [ACC,RT] = GetBehave(history,thresh)
    Rin = squeeze(history(:,1,:));
    Rout = squeeze(history(:,2,:));
    
    in = (Rin > thresh);
    out = (Rout > thresh);
    for j=1:size(history,1)
        t1 = find(in(j,:));
        t2 = find(out(j,:));
        if ~isempty(t1)
            choice(j) = 1;
            RT(j) = t1(1);
        elseif ~isempty(t2)
            choice(j) = 0;
            RT(j) = t2(1);

        else
            choice(j) = 0;
            RT(j) = length(in(1,:));
        end
    end
    
    ACC = mean(choice);
    RT = mean(RT)./1000;
end