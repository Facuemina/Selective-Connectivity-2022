%Test Topology Mod1
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

data_path = Actual_directory(1:end-3); data_path(length(data_path)+1) = {'Data-Selective-Connectivity-2022'};
data_path(length(data_path)+1) = {'Signal Reinforcement'}; data_path(length(data_path)+1) = {'N 5000'};
data_path = join(data_path,Folder_delimiter{fd_choose});
addpath(genpath(data_path{1}))

%%

for i = 1:5
    load(['LoadingCapacity_C_N5000_Epsilon_',num2str(i),'.mat'])
end

LC_DATA = {LoadingCapacity_C_N5000_Epsilon_1,LoadingCapacity_C_N5000_Epsilon_2,LoadingCapacity_C_N5000_Epsilon_3,LoadingCapacity_C_N5000_Epsilon_4,LoadingCapacity_C_N5000_Epsilon_5};

DATA5000_EPC.GENERAL_Weight_Values = cell(1,4);
DATA5000_EPC.GENERAL_Weight_Distr  = cell(1,4);
DATA5000_EPC.GENERAL_Weight_DistrR = cell(1,4);
DATA5000_EPC.Weight_Distr_PC       = cell(1,4);
DATA5000_EPC.Weight_DistrR_PC      = cell(1,4);
DATA5000_EPC.Weight_Distr_PDF      = cell(1,4);
DATA5000_EPC.Weight_DistrR_PDF     = cell(1,4);

c = zeros(1,10);
N = 5000;

DATA5000_EPC.GENERAL_Weight_Values{1} = cell(1,10);
DATA5000_EPC.GENERAL_Weight_Distr{1}  = cell(1,10);
DATA5000_EPC.GENERAL_Weight_DistrR{1} = cell(1,10);
DATA5000_EPC.Weight_Distr_PC{1}       = cell(1,10);
DATA5000_EPC.Weight_DistrR_PC{1}      = cell(1,10);
DATA5000_EPC.Weight_Distr_PDF{1}      = cell(1,10);
DATA5000_EPC.Weight_DistrR_PDF{1}     = cell(1,10);


% for i = 1:10
%     load(['OUT',num2str(i),'.mat'])
%     pmax = max(OUT.pmax);
%     c(i) = LC_DATA{1}.C(i);
%     W    = LC_DATA{1}.Patron(1:pmax,:)'*LC_DATA{1}.Patron(1:pmax,:);
%     
%     Copt = logical(OUT.CMSq);
%     Crand = logical(OUT.CRand);
%     Cwr = logical(1-eye(N));
%     DATA5000_EPC.GENERAL_Weight_Values{1}{i} = unique(W(Cwr));
% 
%     DATA5000_EPC.GENERAL_Weight_Distr{1}{i}  = zeros(1,length(DATA5000_EPC.GENERAL_Weight_DistrR{1}));
%     DATA5000_EPC.GENERAL_Weight_DistrR{1}{i} = zeros(1,length(DATA5000_EPC.GENERAL_Weight_DistrR{1}));
%     
%     for j = 1:length(DATA5000_EPC.GENERAL_Weight_Values{1}{i})
% 
%         DATA5000_EPC.Weight_Distr_PC{1}{i}(j) = length(find(W(Copt) == DATA5000_EPC.GENERAL_Weight_Values{1}{i}(j)))/length(find(W(Cwr) == DATA5000_EPC.GENERAL_Weight_Values{1}{i}(j)));
%         DATA5000_EPC.Weight_DistrR_PC{1}{i}(j)= length(find(W(Crand) == DATA5000_EPC.GENERAL_Weight_Values{1}{i}(j)))/length(find(W(Cwr) == DATA5000_EPC.GENERAL_Weight_Values{1}{i}(j)));
%         
%         DATA5000_EPC.Weight_Distr_PDF{1}{i}(j) = length(find(W(Copt) == DATA5000_EPC.GENERAL_Weight_Values{1}{i}(j)))/(N*LC_DATA{1}.C(i));
%         DATA5000_EPC.Weight_DistrR_PDF{1}{i}(j)= length(find(W(Crand) == DATA5000_EPC.GENERAL_Weight_Values{1}{i}(j)))/(N*LC_DATA{1}.C(i));
%         
%     end
%   
%     clear OUT
%     disp([1,i])
% end

for k = 1:length(LC_DATA)
    DATA5000_EPC.GENERAL_Weight_Values{k} = cell(1,10);
    DATA5000_EPC.GENERAL_Weight_Distr{k}  = cell(1,10);
    DATA5000_EPC.GENERAL_Weight_DistrR{k} = cell(1,10);
    DATA5000_EPC.Weight_Distr_PC{k}       = cell(1,10);
    DATA5000_EPC.Weight_DistrR_PC{k}      = cell(1,10);
    DATA5000_EPC.Weight_Distr_PDF{k}      = cell(1,10);
    DATA5000_EPC.Weight_DistrR_PDF{k}     = cell(1,10);

%     if k == 3
%         problema = [8,9,10];
%     else 
%         problema = 1:10;
%     end
    for i = 1:10
        load(['OUT',num2str(i),'_',num2str(k),'.mat']);
        pmax = max(OUT.pmax);
        c(i) = LC_DATA{k}.C(i);
        W    = OUT.Patron(1:pmax,:)'*OUT.Patron(1:pmax,:);

        Copt = logical(OUT.CMSq);
        Crand = logical(OUT.CRand);
        Cwr = logical(1-eye(N));
        DATA5000_EPC.GENERAL_Weight_Values{k}{i} = unique(W(Cwr));

        DATA5000_EPC.GENERAL_Weight_Distr{k}{i}  = zeros(1,length(DATA5000_EPC.GENERAL_Weight_DistrR{k}));
        DATA5000_EPC.GENERAL_Weight_DistrR{k}{i} = zeros(1,length(DATA5000_EPC.GENERAL_Weight_DistrR{k}));

        for j = 1:length(DATA5000_EPC.GENERAL_Weight_Values{k}{i})

            DATA5000_EPC.Weight_Distr_PC{k}{i}(j) = length(find(W(Copt) == DATA5000_EPC.GENERAL_Weight_Values{k}{i}(j)))/length(find(W(Cwr) == DATA5000_EPC.GENERAL_Weight_Values{k}{i}(j)));
            DATA5000_EPC.Weight_DistrR_PC{k}{i}(j)= length(find(W(Crand) == DATA5000_EPC.GENERAL_Weight_Values{k}{i}(j)))/length(find(W(Cwr) == DATA5000_EPC.GENERAL_Weight_Values{k}{i}(j)));

            DATA5000_EPC.Weight_Distr_PDF{k}{i}(j) = length(find(W(Copt) == DATA5000_EPC.GENERAL_Weight_Values{k}{i}(j)))/(N*LC_DATA{k}.C(i));
            DATA5000_EPC.Weight_DistrR_PDF{k}{i}(j)= length(find(W(Crand) == DATA5000_EPC.GENERAL_Weight_Values{k}{i}(j)))/(N*LC_DATA{k}.C(i));

        end

        clear OUT
        disp([k,i])
    end
end

%% Data Processing
DATA5000_EPC.Mean_Weight_Distr_PC  = cell(1,10);
DATA5000_EPC.Std_Weight_Distr_PC   = cell(1,10);
DATA5000_EPC.Mean_Weight_DistrR_PC = cell(1,10);
DATA5000_EPC.Std_Weight_DistrR_PC  = cell(1,10);
DATA5000_EPC.Mean_Weight_Distr_PDF = cell(1,10);
DATA5000_EPC.Std_Weight_DistrR_PDF = cell(1,10);
DATA5000_EPC.Mean_Weight_Distr_PDF = cell(1,10);
DATA5000_EPC.Std_Weight_DistrR_PDF = cell(1,10);
DATA5000_EPC.Values                = cell(1,10);

Indices_Minimos = zeros(1,10);
Longitud = zeros(5,10);
for i = 1:10 
    for k = 1:length(LC_DATA)
        Longitud(k,i) = length(DATA5000_EPC.GENERAL_Weight_Values{k}{i});
    end
    index = find(round(Longitud(:,i)/2)-Longitud(:,i)/2 ~= 0);
    if isempty(index)
        index = 1:length(LC_DATA);
    end
    [~,Indices_Minimos(i)] = min(Longitud(index,i));
    
%     Vals_Poiss           = DATA5000_EPC.GENERAL_Weight_Values{Indices_Minimos(i)}{i};
    Vals = DATA5000_EPC.GENERAL_Weight_Values{Indices_Minimos(i)}{i}/max(DATA5000_EPC.GENERAL_Weight_Values{Indices_Minimos(i)}{i});
    DiffVals       = diff(Vals);
    Vals_mean      = Vals;
    Vals_mean(1)   = Vals(1) + DiffVals(1)/4;
    Vals_mean(end) = Vals(end) - DiffVals(end)/4;
    index1 = 1:length(LC_DATA); index1 = find(index1 ~= Indices_Minimos(i));
    for j = 2:length(Vals)-1
        Vals_mean(j) = Vals(j)+(DiffVals(j)-DiffVals(j-1))/6;
    end

    DATA5000_EPC.Mean_Weight_Distr_PC{i}  = zeros(1,length(Vals));
    DATA5000_EPC.Std_Weight_Distr_PC{i}   = zeros(1,length(Vals));
    DATA5000_EPC.Mean_Weight_DistrR_PC{i} = zeros(1,length(Vals));
    DATA5000_EPC.Std_Weight_DistrR_PC{i}  = zeros(1,length(Vals));
    DATA5000_EPC.Mean_Weight_Distr_PDF{i} = zeros(1,length(Vals));
    DATA5000_EPC.Std_Weight_DistrR_PDF{i} = zeros(1,length(Vals));
    DATA5000_EPC.Mean_Weight_Distr_PDF{i} = zeros(1,length(Vals));
    DATA5000_EPC.Std_Weight_DistrR_PDF{i} = zeros(1,length(Vals));
        
    for j = 1:length(Vals)
        suma         = DATA5000_EPC.Weight_Distr_PC{Indices_Minimos(i)}{i}(j);
        sumaR        = DATA5000_EPC.Weight_DistrR_PC{Indices_Minimos(i)}{i}(j);
        sumaPDF      = DATA5000_EPC.Weight_Distr_PDF{Indices_Minimos(i)}{i}(j);
        sumaPDF_R    = DATA5000_EPC.Weight_DistrR_PDF{Indices_Minimos(i)}{i}(j);

        for k = 1:length(LC_DATA)
            w     = DATA5000_EPC.GENERAL_Weight_Values{k}{i}/max(DATA5000_EPC.GENERAL_Weight_Values{k}{i});
            
            if j == 1
                index = (w<(Vals(j)+DiffVals(j)/2));
            elseif j == length(Vals)
                index = (w>=(Vals(j)-DiffVals(j-1)/2));
            else
                index = (w>=(Vals(j)-DiffVals(j-1)/2) & w<(Vals(j)+DiffVals(j)/2));
            end
            
            wdist      = DATA5000_EPC.Weight_Distr_PC{k}{i}(index);
            wdistR     = DATA5000_EPC.Weight_DistrR_PC{k}{i}(index);
            wdistPDF   = DATA5000_EPC.Weight_Distr_PDF{k}{i}(index);
            wdistPDF_R = DATA5000_EPC.Weight_DistrR_PDF{k}{i}(index);
            
            suma      = [suma,wdist];
            sumaR     = [sumaR,wdistR];
            sumaPDF   = [sumaPDF,wdistPDF];
            sumaPDF_R = [sumaPDF_R,wdistPDF_R];
        end
        DATA5000_EPC.Mean_Weight_Distr_PC{i}(j)  = mean(suma);
        DATA5000_EPC.Std_Weight_Distr_PC{i}(j)   = std(suma);
        DATA5000_EPC.Mean_Weight_DistrR_PC{i}(j) = mean(sumaR);
        DATA5000_EPC.Std_Weight_DistrR_PC{i}(j)  = std(sumaR);
        DATA5000_EPC.Mean_Weight_Distr_PDF{i}(j) = mean(sumaPDF);
        DATA5000_EPC.Std_Weight_DistrR_PDF{i}(j) = mean(sumaPDF_R);
        DATA5000_EPC.Mean_Weight_Distr_PDF{i}(j) = std(sumaPDF);
        DATA5000_EPC.Std_Weight_DistrR_PDF{i}(j) = std(sumaPDF_R);
        
    end
    
    DATA5000_EPC.Values{i}      = Vals_mean;
%     DATA5000_EPC.Values_dist{i} = Vals_mean_dist;
end
%% PLOTS: DISTRIBUCION P(C=1|W)
figure(1)
sp = 0;
extraInputs = {'fontsize',25};
extraInputs2 = {'fontsize',15};

for i = 1:10%[1,3,5,6,8,9]
  
    sp = sp +1;
%     subplot(2,3,sp)
    figure(sp)
    W_Frac   = DATA5000_EPC.Mean_Weight_Distr_PC{i};
    W_Frac_R = DATA5000_EPC.Mean_Weight_DistrR_PC{i};
    
    Err  = DATA5000_EPC.Std_Weight_Distr_PC{i};
    ErrR = DATA5000_EPC.Std_Weight_DistrR_PC{i};
    
    x = [DATA5000_EPC.Values{i}', fliplr(DATA5000_EPC.Values{i}')];
    xc = W_Frac;
    sc = Err;
    inBetween = [xc+sc, fliplr(xc-sc)];
    h = fill(x, inBetween, 'r');
    set(h,'facealpha',.2,'EdgeColor','r','edgealpha',0.4)
    
    hold on
    
%     x = [DATA5000_EPC.Values{i}', fliplr(DATA5000_EPC.Values{i}')];
%     xc = W_Frac_R;
%     sc = ErrR;
%     inBetween = [xc+sc, fliplr(xc-sc)];
%     h = fill(x, inBetween, 'k');
%     set(h,'facealpha',.2,'EdgeColor','k','edgealpha',0.4)
%         
%     hold on
    
%     plot(DATA5000_EPC.Values{i},W_Frac_R,'k','LineWidth',2)
    plot([-1,1],[1,1]*LoadingCapacity_C_N5000_Epsilon_1.C(i)/5000,'--k','LineWidth',2)
    hold on
    plot(DATA5000_EPC.Values{i},W_Frac,'r','LineWidth',2)
    
    legend({'Optimizada','Aleatoria'},extraInputs2{:});

%     title(['Weight fraction distribution -- c = ',num2str(c(i))])
    xlim([-1,1])
    ylim([0,1])
    
    title(['C = ',num2str(LoadingCapacity_C_N5000_Epsilon_1.C(i)),' - c/N \sim ',num2str(round(LoadingCapacity_C_N5000_Epsilon_1.C(i)/5000,3))])
%     xlabel('W_{ij}/max\{W_{ij}\}')
%     ylabel('frac\{W_{ij}^{opt}\}/frac\{W_{ij}^{rand}\}')
    if i == 8
        xlabel({'W_{k}^{M}'},extraInputs{:})
    end
    
    if i == 1 || i == 6
        ylabel({'P(c_{ij}=1 | W_{ij})'},extraInputs{:})
    end
    
    set(gca,extraInputs{:})
    grid()
end

%% Distribucion J
figure(5)
sp = 0;
extraInputs = {'fontsize',20};
extraInputsL = {'fontsize',15};

for i = 1:10%[1,3,5,6,8,9]
    sp = sp +1;
%     subplot(2,3,sp)
    figure(sp)
    plot(DATA5000_EPC.GENERAL_Weight_Values{1}{i}/max(abs(DATA5000_EPC.GENERAL_Weight_Values{1}{i})),DATA5000_EPC.Weight_Distr_PDF{1}{i},'-r')
    hold on
    
    plot(DATA5000_EPC.GENERAL_Weight_Values{1}{i}/max(abs(DATA5000_EPC.GENERAL_Weight_Values{1}{i})),DATA5000_EPC.Weight_DistrR_PDF{1}{i},'-k')
    hold on

%     plot(Weight_Values{i}/LoadingCapacity_C_N5000_Epsilon_1.C(i),Weight_Distr{i}/(sum(Weight_Distr{i})),'-r')
%     hold on
%     plot(Weight_Values{i}/LoadingCapacity_C_N5000_Epsilon_1.C(i),Weight_DistrR{i}/sum(Weight_DistrR{i}),'-k')
%     hold on
%     plot(Weight_Values{i}/LoadingCapacity_C_N5000_Epsilon_1.C(i),Weight_Distr{i}/(sum(Weight_Distr{i})),'.r')
%     hold on
%     plot(Weight_Values{i}/LoadingCapacity_C_N5000_Epsilon_1.C(i),Weight_DistrR{i}/sum(Weight_DistrR{i}),'.k')
%     hold on

    plot(DATA5000_EPC.GENERAL_Weight_Values{1}{i}/max(abs(DATA5000_EPC.GENERAL_Weight_Values{1}{i})),DATA5000_EPC.Weight_Distr_PDF{1}{i},'.r')
    hold on
    plot(DATA5000_EPC.GENERAL_Weight_Values{1}{i}/max(abs(DATA5000_EPC.GENERAL_Weight_Values{1}{i})),DATA5000_EPC.Weight_DistrR_PDF{1}{i},'-k')
    hold on
    
    title({['C = ',num2str(LoadingCapacity_C_N5000_Epsilon_1.C(i)),' - c/N \sim ',num2str(round(LoadingCapacity_C_N5000_Epsilon_1.C(i)/5000,3))]},extraInputs{:})
%     xlim([-1,1])
%     ylim([0,0.2])
%     xlabel('W_{ij}/max\{W_{ij}\}')
%     ylabel('frac\{W_{ij}^{opt}\}/frac\{W_{ij}^{rand}\}')
    if i == 8
        xlabel({'(2k-p)/c'},extraInputs{:})
%         xlabel({'(2k-p)/c'},extraInputs{:})
    end
    
    if i == 1 || i == 6
        ylabel({'Distribución'},extraInputs{:})
%         ylabel({'Distribución'},extraInputs{:})
    end
    
    if i ~= 1
%         ylim([0,0.025])
        ylim([0,0.025])
    end
%     xlim([-1,1])
    set(gca,extraInputs{:})
    grid()
    legend({'Optimizada','Random'},extraInputsL{:});
end
% 
% 
