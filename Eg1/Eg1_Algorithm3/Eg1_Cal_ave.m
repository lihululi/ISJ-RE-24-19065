function ave_fitness = Eg1_Cal_ave(population, M)
sum_fitness = 0;
for i = 1:M
    sum_fitness = sum_fitness + population.Chromosomes(i).fitness;
end
ave_fitness = sum_fitness/M;