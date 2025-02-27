clc
clear
close all
num_samples = 40;
num_rules = 15;

% Define the function
x = linspace(-5,15,201);
y = zeros(size(x));
y(x>=0 & x<=5) = x(x>=0 & x<=5) - 2;
y(x>5 & x<=10) = 9 - x(x>5 & x<=10);


sample = randsample(length(x) ,num_samples);
sample = sort(sample);
sample_input = x(sample);
sample_output = y(sample);

train = sort(randsample(num_samples,num_samples/2));
test = setdiff(1:num_samples,train);
train_in = sample_input(train);
train_out = sample_output(train);
test_in = sample_input(test);
test_out = sample_output(test);

 
%Q_3_c)

% Partition input universe

input_mf = zeros(num_rules, length(x));

for i = 1 : num_rules
    
    input_mf(i,:) = trimf(x,[train_in(i)-0.5 ,train_in(i), train_in(i)+0.5]);
    figure(2)
    hold on
    plot(x, input_mf)
    title('input membership function')
    xlabel('Input')
    ylabel('Degree of membership')
    xlim auto
    ylim([0 1.05])
end
 hold off


    % Partition output universe
 output_mf = zeros(num_rules, length(y));

 for j = 1 : num_rules
    
    output_mf(j,:) = trimf(x,[train_out(j)-5 ,train_out(j), train_out(j)+5]);
    figure(3)
    hold on
    plot(x, output_mf)
    title('output membership function')
    xlabel('output')
    ylabel('Degree of membership')
    xlim([-3 5])
    ylim([0 1.05])
 end

%% Define the input and output membership functions 

input_fuzzy = zeros(num_rules, length(train_in));
for i = 1 : num_rules
        
    for j = 1: length(train_in)
        
        input_fuzzy(i, j) = train_in(j) * input_mf(i, x == train_in(j));

    end


end

output_fuzzy = zeros(num_rules, length(train_out));
for i = 1 : num_rules
        
    for j = 1: length(train_out)
        
        output_fuzzy(i, j) = train_out(j) * output_mf(i, find(y == train_out(j), 1));

    end


end

coeficient = output_fuzzy ./ input_fuzzy;
coeficient(isnan(coeficient)) = 0;
coeficient(isinf(coeficient)) = 0;


%% applying the test samples


y_estimate = zeros(length(test_in), 1);

for i = 1 : length(test_in)

    y_estimate(i) = sum(test_in(i) .* coeficient(:, i));

end

%% 
figure()
plot(x, y, 'r-', 'LineWidth', 1.5)
hold on
plot(train_in, y_estimate, 'b-', 'LineWidth', 1.5)
hold on
