clear all; close all
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
load('Data_NoiseRed.mat');
Data  = Data_NoiseRed;

%% Plots
extraInputs = {'fontsize',20};
extraInputs2 = {'fontsize',35};
l=0; ll = 2;
%% Enumerate subplots
nIDs = 2;
alphabet = ('A':'Z').';
chars = num2cell(alphabet(1:nIDs));
chars = chars.';
charlbl = strcat('(',chars,')'); % {'(a)','(b)','(c)','(d)'}
%%
colors = [1,0,0;1,0,1;0.5,0,0.5];
h = figure(1);

for k = 1:length(Data)
    subplot(1,2,1)
    errorbar(Data{k}.Conectivity/Data{k}.N,Data{k}.Mean_ImprovRat,Data{k}.Std_ImprovRat,'LineWidth',2,'color',colors(k,:))
    hold on
    
    if k == 2
        subplot(1,2,2)
        errorbar(Data{k}.Conectivity, Data{k}.Mean_LoadingCapacity,Data{k}.Std_LoadingCapacity,'color',colors(k,:),'LineWidth',2)
    end
end

M = 2;
subplot(1,2,1)
hold on
legend({'N = 500','N = 2000','N = 5000'},'Interpreter','latex')


legend boxoff 
set(gca,extraInputs2{:})
xlabel('$c/N$',extraInputs2{:})
ylabel('$p_c^{opt}/p_c^{rand}$',extraInputs2{:})

subplot(1,2,2)
hold on
errorbar(Data{M}.Conectivity,Data{M}.Mean_LoadingCapacity_R,Data{M}.Std_LoadingCapacity_R,'k','LineWidth',2)
set(gca,extraInputs2{:})
xlabel('c',extraInputs2{:})
ylabel('$p_c$',extraInputs2{:})

legend({'Optimized ($\epsilon = 0$)','Random'},'Interpreter','latex')

legend boxoff 
set(groot,'DefaultTextInterpreter','latex');