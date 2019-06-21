function Main()
clear;format short;addpath public;addpath data;addpath adea;
Algorithms = {'ADEA'};
Problems   = {'DTLZ2'};%,'DTLZ2','DTLZ3','DTLZ4'};
Objectives = {3};
RunNum = 1;

for Alg = 1:length(Algorithms)
    for Prob = 1:length(Problems)
        for Obj = 1:length(Objectives)
            Algorithm = Algorithms{Alg};
            Problem   = Problems{Prob};
            M         = Objectives{Obj};
            for Run = 1:RunNum
                MAIN(Algorithm, Problem, M, Run);
                disp([Problem,'_',num2str(M),'_',num2str(Run),' ',datestr(now)]);
            end
        end
    end
end

end