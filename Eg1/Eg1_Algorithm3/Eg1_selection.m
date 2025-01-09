function [parent1, parent2] = Eg1_selection(population)

M = length(population.Chromosomes(:));

if any([population.Chromosomes(:).fitness] < 0 )
    a = 1;
    b = abs( min( [population.Chromosomes(:).fitness] ) );
    Scaled_fitness = a * [population.Chromosomes(:).fitness] + b;
    
    normalized_fitness = [Scaled_fitness] ./ sum([Scaled_fitness]);
else
    normalized_fitness = (1./[population.Chromosomes(:).fitness]) ./ (sum(1./[population.Chromosomes(:).fitness]));
end

[sorted_fintness_values , sorted_idx] = sort(normalized_fitness , 'DESCEND');

for i = 1 : length(population.Chromosomes)
    temp_population.Chromosomes(i).Gene = population.Chromosomes(sorted_idx(i)).Gene;
    temp_population.Chromosomes(i).fitness = population.Chromosomes(sorted_idx(i)).fitness;
    temp_population.Chromosomes(i).normalized_fitness = normalized_fitness(sorted_idx(i));
end

cumsum = zeros(1 , M);

for i = M:-1:1
    for j = i : M
        cumsum(mod(M-i,M)+1) = cumsum(mod(M-i,M)+1) + temp_population.Chromosomes(j).normalized_fitness;
    end
end

R = rand();
parent1_idx = 1;
for i = 1: length(cumsum)
    if R < cumsum(i)
        parent1_idx = mod(M-i,M)+1;
        break;
    end
end

parent2_idx = parent1_idx;

while_loop_stop = 0;

while parent2_idx == parent1_idx
    while_loop_stop = while_loop_stop + 1;
    R = rand();
    if while_loop_stop > 20
        break;
    end
    for i = 1: length(cumsum)
        if R < cumsum(i)
            parent2_idx = mod(M-i,M)+1;
            break;
        end
    end
end

parent1 = temp_population.Chromosomes(parent1_idx);
parent2 = temp_population.Chromosomes(parent2_idx);

end