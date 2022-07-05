clear all; close all;

%% Add Paths
restoredefaultpath

Folder_delimiter{1} = '\'; %Windows
Folder_delimiter{2} = '/'; %Linux
%CHOOSE:
fd_choose = 1;

Actual_directory = split(cd,Folder_delimiter{fd_choose});
fx_path    = Actual_directory(1:end-1); fx_path(length(fx_path)+1) = {'fx'}; fx_path = join(fx_path,Folder_delimiter{fd_choose}); %Add '\fx' folder to path
%optim_path = Actual_directory(1:end-2); optim_path = join(optim_path,Folder_delimiter{fd_choose}); %Add '\Optimization' folder to path
addpath(fx_path{1});% addpath(optim_path{1});

%% Coefficients
% Coefficients estimated from previous runs for a polinomyc function which
% epsilon = 0
% fits P(c/N) = Popt/Prand. It is used for first guess
Coeff = 1.0e+03 *[-1.9450, 6.0410, -4.8990, -3.6990, 9.4120, -7.1130, 2.6780, -0.5289, 0.0523, 0.0026];
Polyn = @(x) polyval(Coeff,x);

%% Define parameters
Dimention   = 100;
Connections = 20;%round(linzpace(20,Dimention-2,30));
e           = 0;%Defines epsilon=0

%% Define Outputs
Outputs = cell(length(Dimention),length(Connections));
Connectivity_Matrices = cell(1,length(Connections));
Pmax    = zeros(length(Dimention),length(Connections));
PRand   = zeros(length(Dimention),length(Connections));
PAT     = cell(size(Dimention));

%% Optimize
for n = 1:length(Dimention)
    
    PAT{n} = sign(rand(Dimention(n),Dimention(n))-0.5); 
    for ind = 1:length(Connections)
        tic
        [PRand(n,ind),C]             = StorageCapacity(Dimention(n),Connections(ind),0,PAT{n}); 
        parameters                   = round([Polyn(Connections(ind)/Dimention(n))*PRand(n,ind),PRand(n,ind)]);
        [Pmax(n,ind),Outputs{n,ind}] = Optimization_SA(Dimention(n),Connections(ind),PAT{n},parameters,e);

        disp(['DONE! - c = ',num2str(Connections(ind)),'-',num2str(ind)])    
        Outputs{n,ind}.CRand = C;
        toc
    end
end

%% Organize Data
Data.Pmax   = Pmax;
Data.PRand  = PRand;
Data.PAT    = PAT;
Data.C      = Connections;
Data.N      = Dimention;

%% Save Data
Folders = {'Data_BothVar','Data_DimentionVar','Data_ConnectivityVar'};

if all([length(Dimention),length(Connections)]>1) || all([length(Dimention),length(Connections)]==1)
    folder_index = 1;
    folder2save = strjoin({cd,'Data',Folders{folder_index}},Folder_delimiter{fd_choose});
elseif length(Dimention)> 1
    folder_index = 2;
    Subfolder = ['C',num2str(Connections(1))];
    folder2save = strjoin({cd,'Data',Folders{folder_index},Subfolder},Folder_delimiter{fd_choose});
elseif length(Connections)> 1    
    folder_index = 3;
    Subfolder = ['N',num2str(Dimention(1))];
    folder2save = strjoin({cd,'Data',Folders{folder_index},Subfolder},Folder_delimiter{fd_choose});
end

if not(isfolder(folder2save))
    mkdir([folder2save,Folder_delimiter{fd_choose},'Connectivity Matrices'])
end

NumData = dir(folder2save);
NumData([NumData.isdir]) = [];
save([folder2save,Folder_delimiter{fd_choose},'Data',num2str(length(NumData)+1),'.mat'],'Data')

for n = 1:length(Dimention)
    for ind = 1:length(Connections)
        OUT = Outputs{n,ind};
        save([folder2save,Folder_delimiter{fd_choose},'Connectivity Matrices',Folder_delimiter{fd_choose},'Data',num2str(length(NumData)+1),'-OUT-N',num2str(Dimention(n)),'-C',num2str(Connections(ind)),'.mat'],'OUT')
    end
end