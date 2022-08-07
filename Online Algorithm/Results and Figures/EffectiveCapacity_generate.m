clear all; close all
%% Add Paths
restoredefaultpath

Folder_delimiter{1} = '\'; %Windows
Folder_delimiter{2} = '/'; %Linux
%CHOOSE:
fd_choose = 1;

Actual_directory = split(cd,Folder_delimiter{fd_choose});
fx_path = Actual_directory(1:end-2); fx_path(length(fx_path)+1) = {'fx'}; fx_path = join(fx_path,Folder_delimiter{fd_choose}); %Add '\fx' folder to path
addpath(fx_path{1});
data_path = Actual_directory(1:end-2); data_path(end+1:end+3) = {'Data';'Online Algorithm';'Effective_Capacity'};
data_path = join(data_path,Folder_delimiter{fd_choose}); addpath(data_path{1});
%%
for i = 1:4
    load([direc,'\LoadingCapacity_C_N2000_Epsilon_',num2str(i),'.mat'])
end
DATA = {LoadingCapacity_C_N2000_Epsilon_1,LoadingCapacity_C_N2000_Epsilon_2,...
    LoadingCapacity_C_N2000_Epsilon_3,LoadingCapacity_C_N2000_Epsilon_4};

PeffData.peff = zeros(length(DATA),length(DATA{1}.C));
PeffData.peff_rand = zeros(length(DATA),length(DATA{1}.C));
for data = 1:length(DATA)
    PeffData.peff(data,:)      = effective_storage(DATA{data},true);
    PeffData.peff_rand(data,:) = effective_storage(DATA{data},false);
    disp(['Listo! ',num2str(data)])
end
    
save([data_path{1},'PeffData.mat'],'PeffData')
    
%%
function [peff] = effective_storage(Data,opt)
    peff = zeros(1,length(Data.Pmax));
    for i = 1:length(Data.Pmax)
        if opt
            C  = Data.OUTS{i}.CMSq;
            pm = Data.Pmax(i);
        else
            C  = Data.OUTS{i}.CRand;
            pm = Data.PRand(i);
        end
        
        frac = 1; Delta = min(2000-pm,round(pm*1.8));
        while frac>0.9
            lb = pm;
            pm = pm + Delta;
            [~, frac] = overlap(Data.Patron,pm,Data.C(i),C,0);            
        end
        ub = pm;
        Delta = -fix((ub-lb)/2);
        while abs(Delta)>=1
            pm = pm + Delta;
            [~, frac] = overlap(Data.Patron,pm,Data.C(i),C,0);  
            if frac>0.9
                lb = pm;Delta = fix((ub-lb)/2);
            else
                ub = pm;Delta = -fix((ub-lb)/2);
            end
        end
        peff(i) = lb;
        disp([i,max(Data.OUTS{i}.pmax),lb])
    end
end