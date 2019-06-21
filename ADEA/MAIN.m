% Decomposition Main File
function MAIN(Algorithm,Problem,M,Run)
format compact;tic;
global Theta W Zmin
% basic settings
[Generations,N,p1,p2] = P_settings(Problem,M); % N为种群数
Evaluations = Generations*N; % max number of fitness evaluations
delta   = 1; % the probability of choosing parents locally
Tm      = ceil(N/10); % the mating neighborhood size
Theta   = 5*ones(N,1);

% weight vector initialization
[N,Ws] = F_weight(p1,p2,M);
for i = 1:N
    Ws(i,:) = Ws(i,:)./norm(Ws(i,:));
end
W = Ws;
Generations = floor(Evaluations/N);

%calculat neighboring angle for angle normalization
cosineW = W*W';
[scosineW, ~] = sort(cosineW, 2, 'descend');
acosW = acos(scosineW(:,2));
refW = rad2deg(acosW);

% perform the mating neighborhood
B     = pdist2(W,W);
[~,B] = sort(B,2);
B     = B(:,1:Tm);

% population initialization
rand('seed', sum(100 * clock));
[Population,Boundary,Coding] = P_objective('init',Problem,M,N);
FunctionValue = P_objective('value',Problem,M,Population);

% extream point
Zmin = min(FunctionValue,[],1);
NonDominated    = P_sort(FunctionValue,'first')==1;
BFunctionValue  = FunctionValue(NonDominated,:);
BPopulation     = Population(NonDominated,:);
[Zmax, ZmaxInd] = max(BFunctionValue,[],1);
ZmaxFV  = BFunctionValue(ZmaxInd,:);
ZmaxPop = BPopulation(ZmaxInd,:);

for Gene = 1 : Generations
    uFunctionValue = (FunctionValue - repmat(Zmin, [size(FunctionValue,1) 1]));
    uFunctionValue = uFunctionValue./repmat(sqrt(sum(uFunctionValue.^2,2)), [1 M]);
    cosine         = sum(uFunctionValue.*W,2);
    Theta          = 0.06*M*(refW + rad2deg(acos(cosine)));
        
    for sub = 1 : N
        % choose the candidate matingpool
        if rand < delta
            P = B(sub,:);
        else
            P = 1:N;
        end
        % binary tounament selection for matingpool
        [MatingPool] = F_mating(sub, P, Population, FunctionValue, Coding, W);
        % generate an offspring
        Offspring    = P_generator(MatingPool,Boundary,Coding,N);
        OffspringFV  = P_objective('value',Problem,M,Offspring);
        % update the ideal point
        Zmin  = min(Zmin,OffspringFV);
        % the translation
        uOffspringFV = (OffspringFV - repmat(Zmin, [size(OffspringFV,1) 1]));
        % the association
        uOffspringFV = uOffspringFV./repmat(sqrt(sum(uOffspringFV.^2,2)), [1 M]);
        cosineOFF    = uOffspringFV*W'; % calculate the cosine values between each solution and each vector
        [~, maxcidx] = max(cosineOFF, [], 2);
        R = maxcidx;
        % update the population
        if F_scalar(FunctionValue(R,:),R) > F_scalar(OffspringFV,R)
            Population(R,:)    = Offspring;
            FunctionValue(R,:) = OffspringFV;
        end
    end
%%
    % update the extream point
    CombineFV      = [FunctionValue;ZmaxFV];%otherOffspringFV];
    CombinePop     = [Population;ZmaxPop];%otherOffspring];
    NonDominated   = P_sort(CombineFV,'first')==1;
    BFunctionValue = CombineFV(NonDominated,:);
    BPopulation    = CombinePop(NonDominated,:);
    [Zmax, ZmaxInd] = max(BFunctionValue,[],1);
    ZmaxFV  = BFunctionValue(ZmaxInd,:);
    ZmaxPop = BPopulation(ZmaxInd,:);

    if(mod(Gene, ceil(Generations*0.1)) == 0) && Gene~=Generations
        W = Ws;
        W = W.*repmat((Zmax - Zmin)*1.0,N,1);
        for i = 1:N
            W(i,:) = W(i,:)./norm(W(i,:));
        end
        B     = pdist2(W,W);
        [~,B] = sort(B,2);
        B     = B(:,1:Tm);
        cosineW = W*W';
        [scosineW, ~] = sort(cosineW, 2, 'descend');
        acosW = acos(scosineW(:,2));
        refW = rad2deg(acosW);
        ZmaxFV  = [];
        ZmaxPop = [];
    end
end

NonDominated  = P_sort(FunctionValue,'first')==1;
Population    = Population(NonDominated,:);
FunctionValue = FunctionValue(NonDominated,:);

TrueValue = P_objective('true',Problem,M,1000);
plot3(TrueValue(:,1), TrueValue(:,2), TrueValue(:,3), '.');
hold on;
plot3(FunctionValue(:,1), FunctionValue(:,2), FunctionValue(:,3), 'ro','MarkerFace', 'r');
xlim([0 max(TrueValue(:,1)) + 0.1]);
ylim([0 max(TrueValue(:,2)) + 0.1]);
zlim([0 max(TrueValue(:,3)) + 0.1]);
xlabel('f_1', 'FontSize', 14);ylabel('f_2', 'FontSize', 14);zlabel('f_3', 'FontSize', 14);
view(135, 30);
hold off;
filename = [Algorithm,'_ ',Problem,'_ ',num2str(M),'_ ',num2str(Run)];
title(filename);

eval(['save Data/',Algorithm,'_',Problem,'_',num2str(M),'_',num2str(Run),' Population FunctionValue'])
end


