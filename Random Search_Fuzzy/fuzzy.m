function [fitness] = fuzzy(num_samples, num_rules, input_mf_width, output_mf_width)

% Define the function
x = linspace(-5,15,201);
y = zeros(size(x));
y(x>=0 & x<=5) = x(x>=0 & x<=5) - 2;
y(x>5 & x<=10) = 9 - x(x>5 & x<=10);
num_samples = round(num_samples);
num_rules = round(num_rules);
if num_samples>60
    num_samples = randi([30 ,55]);
    if rem(num_samples, 2) == 1  % If the number is already even
          num_samples = num_samples + 1;  % Add 1 to make it even
    end
    if num_rules > num_samples/2
        num_rules = randi([4, (num_samples/2)]);
    end
end

if rem(round(num_samples), 2) == 1  % If the number is already even
          num_samples = round(num_samples) + 1;  % Add 1 to make it even
else
    num_samples = round(num_samples);
end


if num_rules<1
    num_rules = 3;
end

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

if input_mf_width<0.1
    input_mf_width = 0.5;
end
    for i = 1 : num_rules
        
        input_mf(i,:) = trimf(x,[train_in(i)-input_mf_width ,train_in(i), train_in(i)+input_mf_width]);
    
    
        
    %     figure(2)
    %     hold on
    %     plot(x, input_mf)
    %     title('input membership function')
    %     xlabel('Input')
    %     ylabel('Degree of membership')
    %     xlim auto
    % %     ylim([0 1.05])
    % end
    %  hold off
    
    
        % Partition output universe
     output_mf = zeros(num_rules, length(y));
    
if output_mf_width<0.1
    output_mf_width = randi([4 7]);
end

     for j = 1 : num_rules
        
        output_mf(j,:) = trimf(x,[train_out(j)-output_mf_width ,train_out(j), train_out(j)+output_mf_width]);
    %     figure(3)
    %     hold on
    %     plot(x, output_mf)
    %     title('output membership function')
    %     xlabel('output')
    %     ylabel('Degree of membership')
    %     xlim([-3 5])
    %     ylim([0 1.05])
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
end
    %matching the size of y_estimate to y
    y_estimate_rescaled = (interp1(1:length(y_estimate), y_estimate, linspace(1, length(y_estimate), length(y))));
    mse = mean((y - y_estimate_rescaled).^2);  %calculating the mean square error
    
    fitness = 1/mse; %defining fitness function
    
    
    
    %% 
%     figure(1)
%     plot(x, y, 'r-', 'LineWidth', 1.5)
%     hold on
%     plot(train_in, y_estimate, 'b-', 'LineWidth', 1.5)
%     hold on
    
    
end

% disp(fitness)