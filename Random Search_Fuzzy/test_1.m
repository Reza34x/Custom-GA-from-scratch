% Initialize variables
num_trials = 10000;
num_wins_stay = 0;
num_wins_switch = 0;

% Run the simulation
for i = 1:num_trials
    % Set up the cups
    cups = [0 0 0];
    prize = randi(3);
    cups(prize) = 1;
    
    % Make a guess
    guess = randi(3);
    
    % Reveal an empty cup
    empty_cups = find(cups == 0);
    empty_cups(empty_cups == guess) = [];
    if isempty(empty_cups)
        % If there is only one empty cup, reveal that one
        empty_cup = guess;
    else
        % If there are two empty cups, choose one at random to reveal
        empty_cup = empty_cups(randi(length(empty_cups)));
    end
    
    % Switch or stay
    switch_cups = setdiff(1:3, [guess empty_cup]);
    if randi(2) == 1
        % Switch
        guess = switch_cups;
    end
    
    % Check if the guess is correct
    if guess == prize
        if length(switch_cups) == 1
            num_wins_switch = num_wins_switch + 1;
        else
            num_wins_stay = num_wins_stay + 1;
        end
    end
end

% Print the results
fprintf("Wins if staying: %d\n", num_wins_stay);
fprintf("Wins if switching: %d\n", num_wins_switch);
fprintf("Probability of winning if staying: %.2f%%\n", 100 * num_wins_stay / num_trials);
fprintf("Probability of winning if switching: %.2f%%\n", 100 * num_wins_switch / num_trials);