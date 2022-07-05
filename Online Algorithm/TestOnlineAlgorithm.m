%% Test Online Algorithm
clear all; close all;
%% Add Paths
restoredefaultpath

Folder_delimiter{1} = '\'; %Windows
Folder_delimiter{2} = '/'; %Linux
%CHOOSE:
fd_choose = 1;

Actual_directory = split(cd,Folder_delimiter{fd_choose});
fx_path    = Actual_directory(1:end-1); fx_path(length(fx_path)+1) = {'fx'}; fx_path = join(fx_path,Folder_delimiter{fd_choose}); %Add '\fx' folder to path
addpath(fx_path{1});
%%
N = 100;
cinit = 20;
Patterns = sign(rand(N,N)-0.5);
e = [1,0,1];
DeltaP = 10;

[C,cc,ccstd,peff,iterations] = OnlineAlgorithm(cinit,50,Patterns,DeltaP,e);

OnlineData.N = N;
OnlineData.cinit = cinit;
OnlineData.Patterns = Patterns;
OnlineData.e = e;
OnlineData.DeltaP = DeltaP;
OnlineData.C = C;
OnlineData.Cmean = cc;
OnlineData.Cstd = ccstd;
OnlineData.peff = peff;
OnlineData.iterations = iterations;

save(['ePC_100Ret_OnlineData_N',num2str(N),'_c',num2str(cinit),'_Delta',num2str(DeltaP),'.mat'],'OnlineData')