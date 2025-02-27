clc
clear
close all

% Defining the function
f = @(x1, x2) 4 * (3 - x1).^2 .* exp(-x1.^2 - (x2 - 1).^2) - 10 * ((x1 / 4) - x1.^4 - 2 * (x2).^5) .* exp(-(x1).^2 - (x2).^2) - 1/2 * exp(-(x1+2).^2 - (x2).^2);

% Defining the parameters
generation = 50;
population = 20;
mutationRate = 0.02;
crossoverRate = 0.5;
select_point1 = [];
select_point2 = [];
point1 = zeros(generation ,population);
point2 = zeros(generation ,population);


% Selection random points from -3 to 3
x = linspace(-3, 3, 100);
y = linspace(-3, 3, 100);
 

select_point1= randi([1, 100] , [1, population]);
select_point2= randi([1, 100] , [1, population]);
point1(1,:) = x(select_point1);
point2(1,:) = y(select_point2); 

% Defining fitness values
fitness = f(point1(1,:), point2(1,:));
average_fitness = zeros(generation+1, 1);
best_fitness = zeros(generation+1, 1);
worst_fitness = zeros(generation+1, 1);
average_fitness(1) = sum(f(point1(1,:), point2(1,:)))/population;
best_fitness(1) = max(f(point1(1,:), point2(1,:)));
worst_fitness(1) = min(f(point1(1,:), point2(1,:)));
fitness_normal = fitness./max(fitness);

% Convert selected points to binary
bin_select_point1 = dec2bin_custom(select_point1);
bin_select_point2 = dec2bin_custom(select_point2);

point1(1,:) = x(select_point1);
point2(2,:) = y(select_point2);
% GA algorithm
for g = 1:generation
    
    % Evaluate fitness
    fitness = f(point1(g,:), point2(g,:));
    average_fitness(g+1) = sum(fitness)/population;
    best_fitness(g+1) = max(fitness);
    worst_fitness(g+1) = min(fitness);
    

    %selection by roulette wheel selection methode
    for i = 1:population/4
        selected(i) = RoulettWheelSelection(fitness_normal);  %selecting by fittnes
        i = i + 1;
    end
% selected = unique(selected);
parent1 = bin_select_point1(selected,:);
parent2 = bin_select_point2(selected,:);
offspring1 = [];
offspring2 = [];

%crossover

for j = 2:length(selected)
    if rand<crossoverRate
        
        [offspring1(j-1,:), offspring1(j,:)] = crossover(parent1(j-1,:) , parent1(j ,:));
        [offspring2(j-1,:), offspring2(j,:)] = crossover(parent2(j-1,:) , parent2(j ,:));
        j = j+1;
    end
end





    
    
    % Mutation
    for i = 1:size(offspring1, 1)
        offspring1(i,:) = mutation(offspring1(i,:), mutationRate);
        offspring2(i,:) = mutation(offspring2(i,:), mutationRate);
    end
    
[~, sortedfitness] = sort(fitness);

min_fitness_select = sortedfitness(1:size(offspring1, 1));


    
% Update decimal selected points
select_point1_off = bin2dec_custom( offspring1);
select_point1_off(select_point1_off > 100) = randi([70 ,97]);
select_point1_off(select_point1_off < 1) = randi([3 ,20]);

select_point2_off = bin2dec_custom(offspring2);
select_point2_off(select_point2_off > 100) = randi([60 ,95]);
select_point2_off(select_point2_off < 1) = randi([3 ,20]);
    if g<generation
        point1(g+1,:) = point1(g,:);
        point2(g+1,:) = point2(g,:);
    
        for count = 1: length(select_point1_off)
            point1(g+1,min_fitness_select(count)) = x(select_point1_off(count));
            point2(g+1,min_fitness_select(count)) = y(select_point2_off(count));
        end
    end
end

% Plotting results
figure;
plot(0:generation, best_fitness, 'r', 'LineWidth', 2);
hold on;
plot(0:generation, average_fitness, 'b', 'LineWidth', 2);
plot(0:generation, worst_fitness, 'g', 'LineWidth', 2);
legend('Best fitness', 'Average fitness', 'Worst fitness');
xlabel('Generation');
ylabel('Fitness');
title('Fitness over generations');
%% 

figure()

    g = 50; %which generation to show
    color = rand(1, 3);
    [X,Y] = meshgrid(x,y);
    Z = f(X,Y);
    contourf(X, Y, Z, 50);
    colorbar;
    title('f(x) Contour Plot');
    xlabel('x');
    ylabel('y');
    hold on
    scatter(point1(g,:), point2(g,:), 'filled', 'MarkerFaceColor', color);
    title('Points');
    xlabel('x');
    ylabel('y');
    xlim([-3 3]);
    ylim([-3 3]);
    hold on



