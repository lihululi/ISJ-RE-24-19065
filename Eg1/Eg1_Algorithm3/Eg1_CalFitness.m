    function [fitness,R] = Eg1_CalFitness(gene, i)
        [flag, gamma, R] = GAEg1_LMI(gene);
        if flag == 1
            fitness = gamma;
        else
            fitness = 1.5;
        end
        disp(['fitness value of ', num2str(i), '-th individual: ', num2str(fitness)]);
        disp(gene');
    end