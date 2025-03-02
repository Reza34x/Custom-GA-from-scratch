function [offspring1, offspring2] = crossover(parent1, parent2)
    % Determine the length of the input parents
    n = length(parent1);
    
    % Choose a random crossover point
    crossover_point = randi(n-1);
    
    % Perform the crossover
    offspring1 = [parent1(1:crossover_point) parent2(crossover_point+1:end)];
    offspring2 = [parent2(1:crossover_point) parent1(crossover_point+1:end)];
end