function [BestChrom, cgcurve, fitness_ave]  = Eg1_GA(M, num_chrom, Miter, Pc, Pm, Er, visuailzation)

XXX = [-1.10884665343369 -0.0670931692638947	-0.00165470850511661	-0.0141732264101360	0.00257046810537039	1.29421062825608e-06	0.00185006554767976	-0.000580728240569499	1.28567004068222e-05	-0.000245331360413909	0.000394889310680300	0.000393811469296874	0.0400000000000000	0.0400000000000000	0.0400000000000000	0.0400000000000000	-0.179142805896377	-1.16713514378544e-08	-0.00956918399377680	-0.372148270745581	-0.603867712815845	-0.00573791924324016	0.0207986895429489	-0.554651986193544	0.00925229477292762	0.0348750687238875];
ub = 1.25*XXX;
lb = .75*XXX;

lb = [lb(1:24),0.2*lb(25),0.2*lb(26)];

chrom_range = [ub;lb];

cgcurve = zeros(1 , Miter);
fitness_ave = zeros(1 , Miter);

g = 1;
disp(['Generation #' , num2str(g)]);
[ population ] = Eg1_Initialization(M, chrom_range, XXX);

for i = 1 : M
    [population.Chromosomes(i).fitness, population.Chromosomes(i).R] = Eg1_CalFitness(population.Chromosomes(i).Gene(:), i);
end

[max_val , indx] = sort([ population.Chromosomes(:).fitness ] , 'ASCEND');

cgcurve(g) = population.Chromosomes(indx(1)).fitness;

fitness_ave(g) = Eg1_Cal_ave(population, M);

disp(['best_fitness: ',num2str(cgcurve(g))]);
disp(['ave_fitness: ',num2str(fitness_ave(g))]);

for g = 2 : Miter
    disp(['Generation #' , num2str(g)]);
   
    for k = 1: 2: M
        [parent1, parent2] = Eg1_selection(population);
        [child1 , child2] = Eg1_crossover(parent1 ,parent2, Pc, 'single');
        [child1] = Eg1_mutation(child1, Pm, chrom_range);
        [child2] = Eg1_mutation(child2, Pm, chrom_range);
        newPopulation.Chromosomes(k).Gene = child1.Gene;
        newPopulation.Chromosomes(k+1).Gene = child2.Gene;
    end
    
    for i = 1 : M
        [newPopulation.Chromosomes(i).fitness, newPopulation.Chromosomes(i).R] = Eg1_CalFitness(newPopulation.Chromosomes(i).Gene(:), i);
    end
    [ newPopulation ] = Eg1_elitism(population, newPopulation, Er, chrom_range);
    
    population = newPopulation;
    cgcurve(g) = population.Chromosomes(1).fitness;
    disp(['best_fitness: ',num2str(cgcurve(g))]);
    fitness_ave(g) = Eg1_Cal_ave(population, M);
    disp(['ave_fitness: ',num2str(fitness_ave(g))]);
end

BestChrom.Gene = population.Chromosomes(1).Gene;
BestChrom.Fitness = population.Chromosomes(1).fitness;
BestChrom.R = population.Chromosomes(1).R;

if visuailzation == 1
    iter = (1:Miter);
    hold on;
    box on;
    plot(iter, cgcurve, 'b--', 'linewidth', 2);
    plot(iter, fitness_ave, 'r', 'linewidth', 2);
    legend('Best fitness value', 'Average fitness value');
    xlabel('Generation number');
    ylabel('Fitness value');
    grid on;
    hold off;
end

end