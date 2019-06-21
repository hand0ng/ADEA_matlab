function Offspring = P_generator(MatingPool,Boundary,Coding,MaxOffspring)
% This function includes the SBX crossover operator and the polynomial mutatoion operator.
% DE and the polynomial mutatoion operator.
    [N,D] = size(MatingPool);
    ProM = 1/D;
    DisM = 20;
    % SBX or DE
    switch Coding
        case 'Real'
            if nargin < 4 || MaxOffspring < 1 || MaxOffspring > N
                MaxOffspring = N;
            end
            ProC = 1;       
            DisC = 30;   	         
            Offspring = zeros(N,D);
            for i = 1 : 2 : N
                beta = zeros(1,D);
                miu  = rand(1,D);
                beta(miu<=0.5) = (2*miu(miu<=0.5)).^(1/(DisC+1));
                beta(miu>0.5)  = (2-2*miu(miu>0.5)).^(-1/(DisC+1));
                beta = beta.*(-1).^randi([0,1],1,D);
                beta(rand(1,D)>ProC) = 1;
                Offspring(i,:)   = (MatingPool(i,:)+MatingPool(i+1,:))/2+beta.*(MatingPool(i,:)-MatingPool(i+1,:))/2;
                Offspring(i+1,:) = (MatingPool(i,:)+MatingPool(i+1,:))/2-beta.*(MatingPool(i,:)-MatingPool(i+1,:))/2;
            end
            Offspring = Offspring(1:MaxOffspring,:);
            if MaxOffspring == 1
                MaxValue = Boundary(1,:);
                MinValue = Boundary(2,:);
            else
                MaxValue = repmat(Boundary(1,:),MaxOffspring,1);
                MinValue = repmat(Boundary(2,:),MaxOffspring,1);
            end
        case 'DE'
            CR = 1;
            F  = 0.5;
            Parent1 = MatingPool(1, :);
            Parent2 = MatingPool(2, :);
            Parent3 = MatingPool(3, :);
            Site = rand(N/3,D) < CR;
            Offspring(Site) = Parent1(Site) + F*(Parent2(Site)-Parent3(Site));
            MaxOffspring = 1;
            MaxValue = Boundary(1,:);
            MinValue = Boundary(2,:);
    end
    % polynomial mutatoion operator
    k    = rand(MaxOffspring,D);
    miu  = rand(MaxOffspring,D);
    Temp = k<=ProM & miu<0.5;
    Offspring(Temp) = Offspring(Temp)+(MaxValue(Temp)-MinValue(Temp)).*((2.*miu(Temp)+(1-2.*miu(Temp)).*(1-(Offspring(Temp)-MinValue(Temp))./(MaxValue(Temp)-MinValue(Temp))).^(DisM+1)).^(1/(DisM+1))-1);
    Temp = k<=ProM & miu>=0.5;
    Offspring(Temp) = Offspring(Temp)+(MaxValue(Temp)-MinValue(Temp)).*(1-(2.*(1-miu(Temp))+2.*(miu(Temp)-0.5).*(1-(MaxValue(Temp)-Offspring(Temp))./(MaxValue(Temp)-MinValue(Temp))).^(DisM+1)).^(1/(DisM+1)));
    Offspring(Offspring>MaxValue) = MaxValue(Offspring>MaxValue);
    Offspring(Offspring<MinValue) = MinValue(Offspring<MinValue);
    if strcmp(Coding, 'Real')
        Offspring = Offspring(1,:);
    end
end