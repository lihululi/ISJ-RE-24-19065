    function [fitness,R] = Eg2_CalFitness(gene, i)
        [flag, gamma, R] = Eg2_LMI(gene);
        if flag == 1
            fitness = gamma;
        else
            fitness = 0.32;
        end
        disp(['fitness value of ', num2str(i), '-th individual: ', num2str(fitness)]);
        disp(gene');
    end