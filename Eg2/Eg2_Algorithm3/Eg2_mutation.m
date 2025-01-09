function [child] = Eg2_mutation(child, Pm, range)
    Gene_no = length(child.Gene);
    
    for k = 1: Gene_no
        R = rand();
        ub = range(1, k);
        lb = range(2, k);
        if R < Pm
            r = rand();
            pp = rand();
            if r < 0.5
                child.Gene(k) = child.Gene(k) - 0.2 * (abs(child.Gene(k)) * (1 - pp/2));
            else
                child.Gene(k) = child.Gene(k) + 0.2 * (abs(child.Gene(k)) * (1 - pp/2));
            end
            child.Gene(k) = Eg2_ifout(child.Gene(k), ub, lb);
        end
    end
end