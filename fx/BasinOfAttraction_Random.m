%% Basin Of Attraction - Simulated Annealing Optimization
function [BOA,C] = BasinOfAttraction_Random(c,C,Pattern,evaluate,Initial_Error_Perc)

    BOA      = zeros(1,length(Initial_Error_Perc)); 
    for perc = 1:length(Initial_Error_Perc)            
        [~, BOA(perc)] = overlap_err(Pattern'*Pattern.*C/c,Pattern,length(Pattern(:,1)),evaluate,Initial_Error_Perc(perc));
        if BOA(perc) == 0
            break
        end
    end
    
end