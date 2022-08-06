%% Test Basin Of Attraction SA
clear all; close all;

%% Add Paths
restoredefaultpath

Folder_delimiter{1} = '\'; %Windows
Folder_delimiter{2} = '/'; %Linux
%CHOOSE:
fd_choose = 1;

Actual_directory = split(cd,Folder_delimiter{fd_choose});
fx_path    = Actual_directory(1:end-2); fx_path(length(fx_path)+1) = {'fx'}; fx_path = join(fx_path,Folder_delimiter{fd_choose}); %Add '\fx' folder to path
addpath(fx_path{1});

%% Define parameters
N        = 2000;
c        = 1500;
step_per = 0.01;
pat_step = 5;
Num2Sto  = c*5.5;
Pat      = sign(rand(N,N)-0.5);
Patterns_To_Evaluate = round(1:pat_step:Num2Sto);

%% Define Outputs
Outputs = cell(1,3);

for boa = 1:length(Outputs)
    Outputs{boa} = zeros(length(Patterns_To_Evaluate),length(0:step_per:0.6));
    if boa == 1
        e = 0;
        Outputs{boa} = BOA_OptimizedNet(N,c,Patterns_To_Evaluate,Pat,step_per,0,e,boa);    
    else
        %Random connectivity matrix
        C = zeros(N,N); v = [ones(1,c),zeros(1,N-c)];
        for i = 1:N
            index      = find(1:N ~= i);
            C(i,index) = v(randperm(N-1));
        end
        %--------------------
        Outputs{boa} = BOA_RandomNet(c,C,Patterns_To_Evaluate,Pat,step_per,0,boa);   
    end
end

BOA_N2000_SA_C1500.Patterns                = Pat;
BOA_N2000_SA_C1500.parameters.N            = N;
BOA_N2000_SA_C1500.parameters.c            = c;
BOA_N2000_SA_C1500.parameters.Num2Sto      = Num2Sto;
BOA_N2000_SA_C1500.parameters.Patterns_To_Evaluate  = Patterns_To_Evaluate;
BOA_N2000_SA_C1500.parameters.pat_step     = pat_step;
BOA_N2000_SA_C1500.parameters.step_per     = step_per;

BOA_N2000_SA_C1500.NonEpsilon_Optimization = Outputs{1};
BOA_N2000_SA_C1500.Epsilon_Optimization    = Outputs{2};
BOA_N2000_SA_C1500.Random_NonOptimization  = Outputs{3};

save(['BOA_N2000_SA_C',num2str(c),'.mat'],'BOA_N2000_SA_C1500')

%% Funtions
function BOA_Data = BOA_OptimizedNet(N,c,Patterns_To_Evaluate,Pat,step_per,evaluate,e,boa)

    BOA_Data = zeros(length(Patterns_To_Evaluate),length(0:step_per:0.6));
    Initial_Error_Perc = 0:step_per:0.6;
    
    pat = 1;
    [BOA_Data(pat,:),~] = BasinOfAttraction_SA(N,c,Pat(1:Patterns_To_Evaluate(1),:),...
        e*[(Patterns_To_Evaluate(pat))/c,1]',evaluate,Initial_Error_Perc);
    disp([num2str(boa),' -- ',num2str(pat),' -- ',num2str(length(Patterns_To_Evaluate))])

    for pat = 2:length(Patterns_To_Evaluate)
        epsilon = e*(Patterns_To_Evaluate(pat))/c;
        
        index_to_use = (BOA_Data(pat-1,:)>0);
        [BOA_Data(pat,index_to_use),~] = BasinOfAttraction_SA(N,c,Pat(1:Patterns_To_Evaluate(pat),:),epsilon,evaluate,Initial_Error_Perc(index_to_use));

        disp([num2str(boa),' -- ',num2str(pat),' -- ',num2str(length(Patterns_To_Evaluate))])
        if BOA_Data(pat,:) == zeros(size(BOA_Data(pat,:)))
            break
        end
    end
end

function BOA_Data = BOA_RandomNet(c,C,Patterns_To_Evaluate,Pat,step_per,evaluate,boa)
    BOA_Data = zeros(length(Patterns_To_Evaluate),length(0:step_per:0.6));
    Initial_Error_Perc = 0:step_per:0.6;
    pat = 1;
    [BOA_Data(pat,:),~] = BasinOfAttraction_Random(c,C,Pat(1:Patterns_To_Evaluate(pat),:),evaluate,Initial_Error_Perc);

%     while BOA_Data(pat,:) ~= ones(1,length(Initial_Error_Perc))
    for pat = 2:length(Patterns_To_Evaluate)
        
        index_to_use = (BOA_Data(pat-1,:)>0);
        [BOA_Data(pat,index_to_use),~] = BasinOfAttraction_Random(c,C,Pat(1:Patterns_To_Evaluate(pat),:),evaluate,Initial_Error_Perc(index_to_use));
        
        disp([num2str(boa),' -- ',num2str(pat),' -- ',num2str(length(Patterns_To_Evaluate))])
        if BOA_Data(pat,:) == zeros(size(BOA_Data(pat,:)))
            break
        end
        
    end
end