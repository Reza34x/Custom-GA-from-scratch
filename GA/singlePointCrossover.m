function offspring = singlePointCrossover(parent1, parent2)
    crossoverPoint = randi(length(parent1)-1)+1; % random crossover point
    offspring = [parent1(1:crossoverPoint), parent2(crossoverPoint+1:end)];
end