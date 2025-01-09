function [BestChrom, cgcurve, fitness_ave]  = Eg2_GA(M, num_chrom, Miter, Pc, Pm, Er, visuailzation)
XXX = [-0.662308204014592	0.982572272962447	0.706219703967137	0.483678310517031	-0.358831927276244	0.00734569281851505	0.00122809665946992	-0.328639295319528	-0.268120864088827	0.00305978566281266	0.00195621192652683	-0.410516515337592	0.100000000000000	0.100000000000000	0.100000000000000	0.100000000000000	0.100000000000000	0.100000000000000	0.100000000000000	0.100000000000000	-10.3303291206101	-10.4598410744258	-0.0366433203891651	0.0504031312980258	-0.0292485219095927	-0.00226731204214527	-13.0137257922116	-13.1203360568281	-11.4352647572985	-11.6085104654116	0.0649808322330891	-0.0428143602402947	-0.0362262557534450	-0.0232122119562211	-13.9317850730245	-14.1731685071712	0.186316418207241	0.0209640549744969];
ub = 1.3.*XXX;
lb = 0.7.*XXX;
chrom_range = [ub;lb];
cgcurve = zeros(1 , Miter);
fitness_ave = zeros(1 , Miter);

g = 1;
disp(['Generation #' , num2str(g)]);
[ population ] = Eg2_Initialization(M, chrom_range, XXX);

for i = 1 : M
    [population.Chromosomes(i).fitness, population.Chromosomes(i).R] = Eg2_CalFitness(population.Chromosomes(i).Gene(:), i);
end

[max_val , indx] = sort([ population.Chromosomes(:).fitness ] , 'ASCEND');
cgcurve(g) = population.Chromosomes(indx(1)).fitness;
fitness_ave(g) = Eg2_Cal_ave(population, M);

disp(['best_fitness: ',num2str(cgcurve(g))]);
disp(['ave_fitness: ',num2str(fitness_ave(g))]);

for g = 2 : Miter
    disp(['Generation #' , num2str(g)]);
   
    for k = 1: 2: M
        [parent1, parent2] = Eg2_selection(population);
        [child1 , child2] = Eg2_crossover(parent1 ,parent2, Pc, 'single');
        [child1] = Eg2_mutation(child1, Pm, chrom_range);
        [child2] = Eg2_mutation(child2, Pm, chrom_range);
        newPopulation.Chromosomes(k).Gene = child1.Gene;
        newPopulation.Chromosomes(k+1).Gene = child2.Gene;
    end
    
    for i = 1 : M
        [newPopulation.Chromosomes(i).fitness, newPopulation.Chromosomes(i).R] = Eg2_CalFitness(newPopulation.Chromosomes(i).Gene(:), i);
    end
    [ newPopulation ] = Eg2_elitism(population, newPopulation, Er, chrom_range);
    
    population = newPopulation;
    cgcurve(g) = population.Chromosomes(1).fitness;
    disp(['best_fitness: ',num2str(cgcurve(g))]);
    fitness_ave(g) = Eg2_Cal_ave(population, M);
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