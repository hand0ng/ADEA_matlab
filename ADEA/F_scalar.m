%Function to calculate scalar value of PBI
function [ScalarValue] = F_scalar(FunctionValue, sub)
    global Theta W Zmin
    [N, ~] = size(FunctionValue);
    normW   = sqrt(sum(W(sub,:).^2,2));
    normP   = sqrt(sum((FunctionValue-repmat(Zmin,N,1)).^2,2));
    CosineP = sum((FunctionValue-repmat(Zmin,N,1)).*repmat(W(sub,:),N,1),2)./normW./normP;
    d1 = normP.*CosineP;
    d2 = normP.*sqrt(1-CosineP.^2);
    ScalarValue = d1 + Theta(sub)*d2;
end