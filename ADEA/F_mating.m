%Function to perform matingpool by binary tounament selection
function [MatingPool] = F_mating(sub, P, Population, FunctionValue, Coding, W)
    switch Coding
        case 'Real'
%            first parent
            P(P == sub) = [];
            MatingPool = Population(sub,:);
%            second parent
            PreP = P(randperm(length(P),2));
            PrePSV = sum(FunctionValue(PreP,:)*W(sub,:)',2);
            [~, P1] = min(PrePSV);
            Parent = PreP(P1);
            MatingPool = [MatingPool; Population(Parent,:)];
        case 'DE'
            % first parent
            P(P == sub) = [];
            MatingPool = Population(sub,:);
            % second parent   
            PreP = P(randperm(length(P),2));
%             PrePSV = F_scalar(FunctionValue(PreP,:), sub);
            PrePSV = sum(FunctionValue(PreP,:),2);
            [~, P1] = min(PrePSV);
            Parent = PreP(P1);
            P(P == PreP(P1)) = [];
            % third parent
            PreP = P(randperm(length(P),2));
%             PrePSV = F_scalar(FunctionValue(PreP,:), sub);
            PrePSV = sum(FunctionValue(PreP,:),2);
            [~, P2] = min(PrePSV);
            Parent = [Parent, PreP(P2)];
            MatingPool = [MatingPool; Population(Parent,:)];
    end
end