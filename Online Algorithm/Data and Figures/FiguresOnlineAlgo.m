% clear all; close all
% % 
% Data = cell(1,3);
% load('C:\Users\Facundo\Desktop\CONICET\Doctorado\Trabajo\Hopfield\Paper 2022\Online Algo\Data Original\Original_OnlineData_N2000_c1999_Delta10.mat')
% Data{1} = OnlineData;
% load('C:\Users\Facundo\Desktop\CONICET\Doctorado\Trabajo\Hopfield\Paper 2022\Online Algo\Data Original\Original_OnlineData_N2000_c1000_Delta10.mat')
% Data{2} = OnlineData;
% % load('C:\Users\Facundo\Desktop\CONICET\Doctorado\Trabajo\Hopfield\Paper 2022\Online Algo\Data Original\Original_OnlineData_N2000_c400_Delta10.mat')
% % Data{end+1} = OnlineData;
% load('C:\Users\Facundo\Desktop\CONICET\Doctorado\Trabajo\Hopfield\Paper 2022\Online Algo\Data Original\Original_OnlineData_N2000_c20_Delta10.mat')
% Data{3} = OnlineData;
% 
% %load('C:\Users\Facundo\Desktop\CONICET\Doctorado\Trabajo\Hopfield\Paper 2022\Loading Capacity\DATA_GRAFICOS_EPC.mat')
% 
% % colors = [0,0,1;1,0,0;0.5,0,0;0.5,0,1];
% colors = [0,0,1;1,0,0;0.7,0.7,0.3];
% %%
% % figure(1)
% % plot(0,0,'color',colors(1,:),'LineWidth',2);hold on;plot(0,0,'color',colors(length(Data),:),'LineWidth',2)
% % for i = [1,length(Data)]
% %     plot_fill_space((1:length(Data{i}.Cmean)),Data{i}.Cmean,Data{i}.Cstd,colors(i,:),2,true)
% %     %xscale('log')
% %     hold on
% % end
% % xlim([0,12000])
% % ylabel('$\langle c\rangle$','Interpreter','latex')
% % xlabel('Iterations','Interpreter','latex')
% % legend({'$c_0=1999$','$c_0=20$'},'Interpreter','latex')
% % set(gca,'FontSize',25)
% %%
% load('LoadingCapacity_C_N2000_Epsilon_Paper.mat')
% for i = 2:4
%     load(['LoadingCapacity_C_N2000_Epsilon_Paper_',num2str(i)])
% end
% Peff = mean([LoadingCapacity_C_N2000_Epsilon_Paper.Pmax;LoadingCapacity_C_N2000_Epsilon_Paper_2.Pmax;...
%     LoadingCapacity_C_N2000_Epsilon_Paper_3.Pmax;LoadingCapacity_C_N2000_Epsilon_Paper_4.Pmax]);
% Peff_std = std([LoadingCapacity_C_N2000_Epsilon_Paper.Pmax;LoadingCapacity_C_N2000_Epsilon_Paper_2.Pmax;...
%     LoadingCapacity_C_N2000_Epsilon_Paper_3.Pmax;LoadingCapacity_C_N2000_Epsilon_Paper_4.Pmax]);
figure(1)
% for i = 1:length(Data)
%     plot(0,0,'color',colors(i,:),'LineWidth',2)
%     hold on
% end
% errorbar(DATA_GRAFICOS_EPC{2}.Conectivity,DATA_GRAFICOS_EPC{2}.Mean_LoadingCapacity,DATA_GRAFICOS_EPC{2}.Std_LoadingCapacity,'color',[0.85,0.325,0.098],'LineWidth',2)
% hold on
% errorbar(DATA_GRAFICOS_EPC{2}.Conectivity,DATA_GRAFICOS_EPC{2}.Mean_LoadingCapacity_R,DATA_GRAFICOS_EPC{2}.Std_LoadingCapacity_R,'color','k','LineWidth',2)
IND = unique([1:2:length(DATA_GRAFICOS_EPC{2}.Conectivity),30]);
errorbar(LoadingCapacity_C_N2000_Epsilon_Paper.C,Peff,Peff_std,'color',[0.85,0.325,0.098],'LineWidth',2)
hold on
errorbar(DATA_GRAFICOS_EPC{2}.Conectivity(IND),mean(PeffData.peff_rand(:,IND)),std(PeffData.peff_rand(:,IND)),'color','k','LineWidth',2)
% errorbar(LoadingCapacity_C_N2000_Epsilon_Paper_2.C,peff2,std_peff2)
% hold on
% errorbar(DATA{1}.C(IND),mean(PeffData.peff_rand(:,IND)),std(PeffData.peff_rand(:,IND)),'color','k','LineWidth',2)
for i = 1:7
    plot(0,0,'color',[1,1,1])
    hold on
end
E = 0;
for i = 1:length(Data)
    plot(Data{i}.Cmean(1:end-E),Data{i}.peff(1:end-E),'color',colors(i,:),'LineWidth',3)
%     [~,iM] = max(Data{i}.peff);iM=iM-round(0.5*iM);
%     plot(Data{i}.Cmean(1:iM),Data{i}.peff(1:iM),'color',colors(i,:),'LineWidth',2)
    hold on
end
% errorbar(LoadingCapacity_C_N2000_Epsilon_Paper.C,Peff,Peff_std,'color',[0.85,0.325,0.098],'LineWidth',2)
% hold on
% errorbar(DATA_GRAFICOS_EPC{2}.Conectivity(IND),mean(PeffData.peff_rand(:,IND)),std(PeffData.peff_rand(:,IND)),'color','k','LineWidth',2)
% legend({'$c_0 = 1999$','$c_0 = 1000$','$c_0 = 20$','Optimized ($\epsilon=\alpha$)','Random'},'Interpreter','latex','NumColumns',2)
legend({'Optimized ($\epsilon=\alpha$)','Random','','','','','','','','$c_0 = 1999$','$c_0 = 1000$','$c_0 = 20$'},'Interpreter','latex','NumColumns',4)
% legend({'Optimized ($\epsilon=\alpha$)','Random','','$c_0 = 1999$','$c_0 = 1000$','$c_0 = 20$'},'Interpreter','latex','NumColumns',2)
% legend({'$c_0 = 1999$','$c_0 = 1000$','$c_0 = 20$','Optimized ($\epsilon=\alpha$)','Random'},'Interpreter','latex','NumColumns',2)
xlabel('$\langle c\rangle$','Interpreter','latex')
ylabel('$p_{eff}$','Interpreter','latex')
ylim([0,1100])
legend boxoff
set(gca,'FontSize',25)


function plot_fill_space(X_Data,Y_Data,Y_err_Data,color,LW,log_true)
        %Fill space between errorbars
        x = [X_Data, fliplr(X_Data)];
        inBetween = [Y_Data+Y_err_Data, fliplr(Y_Data-Y_err_Data)];
        h_f = fill(x, inBetween, color);
        set(h_f,'facealpha',.2,'EdgeColor',color,'edgealpha',0.4)
        hold on
        plot(X_Data,Y_Data,'color',color,'linewidth',LW)
end