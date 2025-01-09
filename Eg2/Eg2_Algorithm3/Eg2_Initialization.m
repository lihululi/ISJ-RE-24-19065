function population = Eg2_Initialization(M, chrom_range, XXX0)
s = rng(1019);

for i = 1 : M
    for j = 1 : length(XXX0)
        population.Chromosomes(i).Gene(j) = [rand()*(chrom_range(1, j)-chrom_range(2, j))+chrom_range(2, j)];
    end
end