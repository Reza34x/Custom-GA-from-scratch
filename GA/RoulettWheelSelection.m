function i = RoulettWheelSelection(p)
    r = rand * sum(p); % Find a random number between [0, max(cumsum)=sum]
    c = cumsum(p); % Calculate the cumulative summation of each fitness
    i = find(r <=c, 1, 'first'); % Find the index than which the random number is lower
end



