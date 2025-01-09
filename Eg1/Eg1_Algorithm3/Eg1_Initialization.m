function population = Eg1_Initialization(M, chrom_range, XXX0)
Step = M+1;
s = rng(1202);

for i = 1 : M
    if mod(i,Step)==0
        for j = 1 : length(XXX0)
            population.Chromosomes(i).Gene(j) = XXX0(j);
        end
    else
        for j = 1 : length(XXX0)
            population.Chromosomes(i).Gene(j) = [rand()*(chrom_range(1, j)-chrom_range(2, j))+chrom_range(2, j)];
        end
    end
end