function [new_solution] = simulated_annealing(solution, fitness, f, x, y)
    temperature = 1.0;
    alpha = 0.99;
    k = 1;
    new_solution = solution;
    new_fitness = fitness;
    while temperature > 1e-8
        % Generate a random neighbor solution
        neighbor = [x(randi(length(x), 1)), y(randi(length(y), 1))];
        new_fitness = f(neighbor(1), neighbor(2));
        % Compute the acceptance probability
        delta = new_fitness - fitness;
        if delta < 0 || exp(-delta/temperature) > rand()
            new_solution = neighbor;
            fitness = new_fitness;
        end
        % Reduce the temperature
        temperature = alpha * temperature;
        k = k + 1;
    end
end