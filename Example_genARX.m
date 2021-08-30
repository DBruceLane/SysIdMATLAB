close all;
clc;
clear;

T = importdata('ARX.csv');
u = T(:,1); % Entrada
y = T(:,2); % Saida

% Numero de regressores de entrada e saida
ryMax = 10;
ruMax = 10;

bestFit = -Inf;
i = 1;
arrayFit = zeros(ryMax+ruMax,1);
for ry = 1:ryMax
    for ru = 1:ruMax
        [fitness,MSE,y_hat,Psi,Theta,ksi] = genARX(y,u,ry,ru,0);
        arrayFit(i) = fitness;
        i = i+1;
        
        if fitness>bestFit
            bestFit = fitness;
            bestEst = y_hat;
            bestRy = ry;
            bestRu = ru;
        end
    end
end

% Fitness vs regressores
figure;

yyaxis left;
plot(arrayFit,'LineWidth',2);
hold on;
plot((bestRy-1)*ruMax+bestRu,bestFit,'k*','LineWidth',2);
title('Fitness vs numero de regressores');

yyaxis right;
plotRy = ones(1,ruMax);
for i = 2:ruMax
    plotRy = [plotRy i*ones(1,ruMax)];
end
plotRu = 1:1:ryMax;
for i = 2:ruMax
    plotRu = [plotRu 1:1:ryMax];
end

plot(plotRy,'--r');
hold on;
plot(plotRu,'--g');
ylabel('Ry e Ru');
legend('Fitness','Maior fitness','Regressores y','Regressores u');

% Melhor Estimativa
figure;
plot(y);
hold on;
plot(bestEst);
legend('Saida Real','Saida Estimada');
title(['Melhor Estimativa ARX Simulada' newline num2str(bestRy) 'y ' num2str(bestRu) 'u / Fitness=' num2str(bestFit)]);



