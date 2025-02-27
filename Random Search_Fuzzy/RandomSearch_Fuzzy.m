clc
clear

% Define the range of parameters.

num_samples     = 10:4:60;                %    must be even
num_rules       = 2:max(num_samples)/2;   %    must be <= num_sample/2 
input_mf_width  = linspace(0.1, 2, 10);
output_mf_width = linspace(2, 7, 10);
variables = zeros(length(num_samples),length(num_rules),length(input_mf_width),length(output_mf_width));

best = [];
%% Algorithm
finalvalue = [];
for count = 1:80

% Select Bias
b = 3;

% Initial Conditions (where do we start?)
pts(1,1) = num_samples(randi(numel(num_samples)));        % initial x position
pts(1,2) = num_rules(randi(numel(num_rules)));            % initial y position
while pts(1,2)>pts(1,1)/2                                   
      pts(1,2) = num_rules(randi(numel(num_rules)));
    if pts(1,2)<=pts(1,1)/2
        continue
    end
end
pts(1,3) = input_mf_width(randi(numel(input_mf_width)));  % initial z position
pts(1,4) = output_mf_width(randi(numel(output_mf_width)));% initial w position 
pts(1,5) = fuzzy(pts(1,1),pts(1,2),pts(1,3),pts(1,4));                    
%% 

% Set the convergence criteria
epsilon = 0.01;

% Algorithm
err_x = 5;
err_y = 5;


for i = 1:20 
    % Select random deltas for x and y
    dx = randi(10) * 2;
    dy = randi(2);
    dz = rand(1);
    dw = 4*rand(1);
    if b==0
        b =1;
    end
    % Check positive x direction
    if  pts(i,2)<(pts(i,1)+b*dx)/2 && pts(i,1)+b*dx>0 && fuzzy(pts(i,1)+b*dx,pts(i,2), pts(i,3),pts(i,4)) > fuzzy(pts(i,1),pts(i,2), pts(i,3),pts(i,4))
        pts(i+1,1) = pts(i,1) + b*dx;
        pts(i+1,2) = pts(i,2);
        pts(i+1,3) = pts(i,3);
        pts(i+1,4) = pts(i,4);
        b =0.2 * b + 0.4 * dx;
    % Check negative x direction
    elseif pts(i,2)<(pts(i,1)-b*dx)/2 && pts(i,1)-b*dx>0 && fuzzy(pts(i,1)-b*dx,pts(i,2), pts(i,3),pts(i,4)) > fuzzy(pts(i,1),pts(i,2), pts(i,3),pts(i,4))  
        pts(i+1,1) = pts(i,1) - b*dx;
        pts(i+1,2) = pts(i,2);
        pts(i+1,3) = pts(i,3);
        pts(i+1,4) = pts(i,4);
        b = 0.2 * b - 0.4 * dx;
     % Check positive y direction
    elseif  pts(i,2)+b*dx<(pts(i,1))/2 && pts(i,2)+b*dy>0 && fuzzy(pts(i,1),round(pts(i,2)+b*dy), pts(i,3),pts(i,4)) > fuzzy(pts(i,1),pts(i,2), pts(i,3),pts(i,4))
        pts(i+1,1) = pts(i,1);
        pts(i+1,2) = pts(i,2) + b*dy;
        pts(i+1,3) = pts(i,3);
        pts(i+1,4) = pts(i,4);
        b = 0.2 * b + 0.4 * dy;
    % Check negative y direction
    elseif pts(i,2)-b*dx<(pts(i,1))/2 && pts(i,2)-b*dy>0 && fuzzy(pts(i,1),round(pts(i,2)-b*dy), pts(i,3),pts(i,4)) > fuzzy(pts(i,1),pts(i,2), pts(i,3),pts(i,4))
        pts(i+1,1) = pts(i,1);
        pts(i+1,2) = pts(i,2) - b*dy;
        pts(i+1,3) = pts(i,3);
        pts(i+1,4) = pts(i,4);
        b = 0.2 * b - 0.4 * dy;
         % Check positive z direction
    elseif pts(i,3)+b*dz>0 && fuzzy(pts(i,1),pts(i,2), pts(i,3)+b*dz,pts(i,4)) > fuzzy(pts(i,1),pts(i,2), pts(i,3),pts(i,4))
        pts(i+1,1) = pts(i,1);
        pts(i+1,2) = pts(i,2);
        pts(i+1,3) = pts(i,3) + b*dz;
        pts(i+1,4) = pts(i,4);
        b = 0.2 * b + 0.4 * dz;
    % Check negative z direction
    elseif pts(i,3)-b*dz>0 && fuzzy(pts(i,1),pts(i,2), pts(i,3)-b*dz,pts(i,4)) > fuzzy(pts(i,1),pts(i,2), pts(i,3),pts(i,4))
        pts(i+1,1) = pts(i,1);
        pts(i+1,2) = pts(i,2);
        pts(i+1,3) = pts(i,3) - b*dz;
        pts(i+1,4) = pts(i,4);
        b = 0.2 * b - 0.4 * dz;
         % Check positive w direction
    elseif pts(i,4)+b*dz>0 && fuzzy(pts(i,1),pts(i,2), pts(i,3),pts(i,4)+b*dw) > fuzzy(pts(i,1),pts(i,2), pts(i,3),pts(i,4))
        pts(i+1,1) = pts(i,1);
        pts(i+1,2) = pts(i,2);
        pts(i+1,3) = pts(i,3);
        pts(i+1,4) = pts(i,4) + b*dw;
        b = 0.2 * b + 0.4 * dw;
    % Check negative w direction
    elseif pts(i,4)-b*dz>0 && fuzzy(pts(i,1),pts(i,2), pts(i,3),pts(i,4)-b*dw) > fuzzy(pts(i,1),pts(i,2), pts(i,3),pts(i,4))
        pts(i+1,1) = pts(i,1);
        pts(i+1,2) = pts(i,2);
        pts(i+1,3) = pts(i,3) ;
        pts(i+1,4) = pts(i,4)- b*dw;
        b = 0.2 * b - 0.4 * dw;
    
    else % No improvement, reduce bias
        b = round(0.5 * b);
        pts(i+1,1) = pts(i,1);
        pts(i+1,2) = pts(i,2);
        pts(i+1,3) = pts(i,3);
        pts(i+1,4) = pts(i,4);
    end

        %check if parameters are outside their bouneries
    if  abs(pts(i,1))> 190
        pts(i,1) = 2 * randi([80 90]);
    elseif  abs(pts(i,1))< 2
        pts(i,1) = 2 * randi([3 10]);
    end
      
%     elseif  abs(pts(i,2))> 190
%         pts(i,2) = 2 * randi([80 90]);
%     elseif  abs(pts(i,2))> 190
%         pts(i,2) = 2 * randi([80 90]);
    if  abs(pts(i,3))> 4
        pts(i,3) = 2 * rand();
    elseif  abs(pts(i,3))< 0
        pts(i,3) =  rand();
    end
    if  abs(pts(i,4))> 8
        pts(i,4) = randi([5 8]);
    elseif  abs(pts(i,4))< 2
        pts(i,4) =randi([3 7]);
    end

    % Break loop if going too long or outside of the area
    if pts(i, 5) < 0.1 && i==2
        break;
    
    end
     
    

    % Update the point
    pts(i+1,5) = fuzzy(pts(i+1,1),pts(i+1,2),pts(i+1,3),pts(i+1,4));

    % Update difference between consecutive points
    err_x = abs(pts(i+1,1) - pts(i,1));
    err_y = abs(pts(i+1,2) - pts(i,2));
    
    % Increment i (move to next row)
    disp("pts is")
   pts(i,:)
   disp("i is")
    disp(i)
    disp("count is")
    disp(count)
    
    
% end

max_fitness_select = find(pts(:, 5) == max(pts(:, 5)), 1);
best(count, :) = pts(max_fitness_select, :);


end

% finalvalue(1:count)
% disp('the maximum point is:')
% disp(max(finalvalue))
% disp('the average of answers is:')
% disp(sum(finalvalue)/count)
disp(pts(i,5))

end










