function i = RoulettWheelSelection(p)
    r = rand * sum(p);
    c = cumsum(p);
    i = find(r <=c, 1, 'first');
end



