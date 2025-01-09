function [ newPopulation2 ] = Eg2_elitism(population , newPopulation, Er, chrom_range)

M = length(population.Chromosomes);
Elite_no = round(M * Er);

[max_val , indx] = sort([ population.Chromosomes(:).fitness ] , 'ASCEND');

[max_valn , indxn] = sort([ newPopulation.Chromosomes(:).fitness ] , 'ASCEND');

for k = 1 : Elite_no
    bestPopulation.Chromosomes(2*k-1).Gene = population.Chromosomes(indx(k)).Gene;
    bestPopulation.Chromosomes(2*k-1).fitness = population.Chromosomes(indx(k)).fitness;
    bestPopulation.Chromosomes(2*k-1).R = population.Chromosomes(indx(k)).R;
    bestPopulation.Chromosomes(2*k).Gene = newPopulation.Chromosomes(indxn(k)).Gene;
    bestPopulation.Chromosomes(2*k).fitness = newPopulation.Chromosomes(indxn(k)).fitness;
    bestPopulation.Chromosomes(2*k).R = newPopulation.Chromosomes(indxn(k)).R;
end
[max_valb , indxb] = sort([ bestPopulation.Chromosomes(:).fitness ] , 'ASCEND');
for k = 1 : Elite_no
    newPopulation2.Chromosomes(k).Gene = bestPopulation.Chromosomes(indxb(k)).Gene;
    newPopulation2.Chromosomes(k).fitness  = bestPopulation.Chromosomes(indxb(k)).fitness;
    newPopulation2.Chromosomes(k).R  = bestPopulation.Chromosomes(indxb(k)).R;
end

for k = Elite_no + 1 :  M
    if k <= round(M*(1-Er))
        newPopulation2.Chromosomes(k).Gene = newPopulation.Chromosomes(indxn(k)).Gene;
        newPopulation2.Chromosomes(k).fitness = newPopulation.Chromosomes(indxn(k)).fitness;
        newPopulation2.Chromosomes(k).R = newPopulation.Chromosomes(indxn(k)).R;
    else
        for j = 1 : length(population.Chromosomes(indx(k)).Gene(:))
            newPopulation2.Chromosomes(k).Gene(j) = [rand()*(chrom_range(1, j)-chrom_range(2, j))+chrom_range(2, j)];
        end
        [newPopulation2.Chromosomes(k).fitness, newPopulation2.Chromosomes(k).R] = Eg2_CalFitness(newPopulation2.Chromosomes(k).Gene(:), k);
    end
end