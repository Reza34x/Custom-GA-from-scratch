function y = mutation(chromosome, mu)
% chromosome - a binary vector representing the chromosome
% mutationRate - the probability of each bit in the chromosome to be mutated

    y = chromosome; % create a copy of the input chromosome
    
    for i = 1:length(chromosome)
        if rand < mu % check if the bit should be mutated
            y(i) = ~chromosome(i); % flip the bit
        end
    end
    
end