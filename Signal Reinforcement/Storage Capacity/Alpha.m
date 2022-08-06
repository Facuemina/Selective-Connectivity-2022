%% ALPHA
clear all; close all;
%% Add Paths
restoredefaultpath

Folder_delimiter{1} = '\'; %Windows
Folder_delimiter{2} = '/'; %Linux
%CHOOSE:
fd_choose = 1;

Actual_directory = split(cd,Folder_delimiter{fd_choose});
data_path = Actual_directory(1:end-2); data_path(length(data_path)+1) = {'Data'}; data_path(length(data_path)+1) = {'Storage Capacity'}; 
data_path = join(data_path,Folder_delimiter{fd_choose}); %Add '\fx' folder to path
addpath(data_path{1});
%% Load Data
load('Data_SignalReinf.mat');
load('Data_NoiseRed.mat');
%% Original
errorbar(0,0,0,'color',[0.85,0.325,0.098],'LineWidth',2);hold on
errorbar(0,0,0,'m','LineWidth',2);errorbar(0,0,0,'k','LineWidth',2);
errorbar(Data_SignalReinf{2}.Conectivity/2000,...
    Data_SignalReinf{2}.Mean_LoadingCapacity./Data_SignalReinf{2}.Conectivity,...
    Data_SignalReinf{2}.Std_LoadingCapacity./Data_SignalReinf{2}.Conectivity,...
    'color',[0.85,0.325,0.098],'linewidth',2)
hold on
errorbar(Data_NoiseRed{2}.Conectivity/2000,...
    Data_NoiseRed{2}.Mean_LoadingCapacity./Data_NoiseRed{2}.Conectivity,...
    Data_NoiseRed{2}.Std_LoadingCapacity./Data_NoiseRed{2}.Conectivity,...
    'm','linewidth',2)


errorbar(Data_SignalReinf{2}.Conectivity/2000,...
    Data_SignalReinf{2}.Mean_LoadingCapacity_R./Data_SignalReinf{2}.Conectivity,...
    Data_SignalReinf{2}.Std_LoadingCapacity_R./Data_SignalReinf{2}.Conectivity,...
    'k','linewidth',2)
    
ylabel('$\alpha_c$','Interpreter','latex')
xlabel('$c/N$','Interpreter','latex')
lgd=legend({'Optimization ($\epsilon=p$)','Optimization ($\epsilon=0$)','Random'},'Interpreter','latex');
legend boxoff
set(gca,'fontsize',25)