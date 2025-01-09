function [child1 , child2] = Eg1_crossover(parent1 , parent2, Pc, crossoverName)

switch crossoverName
    case 'single'
        Gene_no = length(parent1.Gene);
        ub = Gene_no - 1;
        lb = 1;
        Cross_P = round ( (ub - lb) *rand() + lb );
        Cross_L = round((ub - lb) * rand() * 0.15 + lb);
        if (Cross_P + Cross_L) > Gene_no
            Cross_L = Gene_no - Cross_P;
        end
        temp1 = parent1.Gene;
        temp2 = parent2.Gene;
        for i = Cross_P : (Cross_P + Cross_L)
            pp = rand();
            if pp > 0.5
                pp = pp * 0.3;
            else
                pp = pp * 0.6;
            end
            temp1(i) = parent1.Gene(i) * (1 - pp) + parent2.Gene(i) * pp;
            temp2(i) = parent2.Gene(i) * (1 - pp) + parent1.Gene(i) * pp;
        end

        child1.Gene = temp1;
    
        child2.Gene = temp2;
        
    case 'double'
        Gene_no = length(parent1);

        ub = length(parent1.Gene) - 1;
        lb = 1;
        Cross_P1 = round (  (ub - lb) *rand() + lb  );
        
        Cross_P2 = Cross_P1;
        
        while Cross_P2 == Cross_P1
            Cross_P2 = round (  (ub - lb) *rand() + lb  );
        end
        
        if Cross_P1 > Cross_P2
            temp =  Cross_P1;
            Cross_P1 =  Cross_P2;
            Cross_P2 = temp;
        end

        Part1 = parent1.Gene(1:Cross_P1);
        Part2 = parent2.Gene(Cross_P1 + 1 :Cross_P2);
        Part3 = parent1.Gene(Cross_P2+1:end);
        
        child1.Gene = [Part1 , Part2 , Part3];
        
        
        Part1 = parent2.Gene(1:Cross_P1);
        Part2 = parent1.Gene(Cross_P1 + 1 :Cross_P2);
        Part3 = parent2.Gene(Cross_P2+1:end);
        
        child2.Gene = [Part1 , Part2 , Part3];
end

R1 = rand();

if R1 <= Pc
    child1 = child1;
else
    child1 = parent1;
end

R2 = rand();

if R2 <= Pc
    child2 = child2;
else
    child2 = parent2;
end

end