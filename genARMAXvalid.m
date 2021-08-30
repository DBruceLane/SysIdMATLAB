function [MSE,y_hat,Psi_EMQ,Theta_EMQ,ksi2] = genARMAXvalid(yOriginal,uOriginal,ry,ru,rn,pEst,mostrarEstimativa)
% Generalizacao de modelagem ARMAX com validação
nOriginal = length(uOriginal);

y = yOriginal(1:pEst*nOriginal);
u = uOriginal(1:pEst*nOriginal);
% pEst (0-1) é a porcentagem dos dados que serão usados pra estimar
% Os dados para validar serão 1-percenEst

n = length(u);

%% Validar Resíduo
var_ksi = 100; %Só para inicializar o loop while
iteracoes = 0; %Contador de iterações para convergir o ksi
while var_ksi>0.5     
    %Se o menor índice na equação de y(k) é (k-1), o menor k possível é 2 (pro Matlab iniciar índice em 1)
    %Se fosse o caso kmin = 2
    if ru>= rn && ru>= ry
        kmin = 1+ru;
    elseif ry>=rn && ry>=ru
        kmin = 1+ry;
    else
        kmin = 1+rn;
    end
    
    Y = y(kmin:n);
    
    if iteracoes==0
        [~,~,Psi,Theta] = genARX2(y,u,ry,ru,0);
        ksi = Y - Psi*Theta; %Calculo do primeiro residuo pelo modelo ARX  
    end

    Psiy = zeros(1+n-kmin,ry); %Elementos de saída na matriz Psi
    for i=1:ry
            Psiy(1:1+n-kmin,i) = y(kmin-i:n-i); 
    end

    Psiu = zeros(1+n-kmin,ru); %Elementos de entrada na matriz Psi
    for i=1:ru
           Psiu(1:1+n-kmin,i) = u(kmin-i:n-i);     
    end
    
    Psir = ksi;

    Psi_EMQ = [Psiy Psiu Psir];
    Theta_EMQ = inv(Psi_EMQ'*Psi_EMQ)*Psi_EMQ'*Y; %Pseudo Inversa
    
    y_hat = Psi_EMQ*Theta_EMQ;
    ksi = Y - y_hat;
    var_ksi = var(ksi);
    iteracoes = iteracoes+1;
    if iteracoes>200
        break
    end
end

%% Free Run Simulation com Resíduo Validado

% Para Validação da Simulação é usado o restante dos dados
if pEst ~= 1
    y = yOriginal(pEst*nOriginal:end);
    u = uOriginal(pEst*nOriginal:end);
else
    y = yOriginal;
    u = uOriginal;
end
n = length(u);

if ru>=(ry+2)
    kmin = ru;  
    %Quando ha 2 regressores a mais de entrada que saida, ex: ry = 1 e ru = 3
    %y(k) = y(k-1) + u(k) + u(k-1) + u(k-2) 
    %Nesse caso observa-se o ultimo elemento de entrada, kmin seria 3
else
    kmin = 1+ry;
    %Nos outros casos observa-se o ultimo elemento de saida
end


clear y_hat
y_hat = zeros(1,ry);
for i=1:ry
    y_hat(i) = y(i); %Condições de Contorno
end

parcial = zeros(1,ry+ru); %Parciais de y_hat
%Exemplo de equação pra mostrar os índices do for
%y_hat(k) = Theta(j=1)*y(k-1) + Theta(j=2)*y(k-2) + Theta(j=3)*u(k) + Theta(j=4)*u(k-1)
for k=kmin:n  
   y_hat(k) = 0;
   for j=1:ry+ru
       if j>ry+ru
           parcial(j) = Theta_EMQ(j)*ksi(k-j+ry); %Parcial dependente do residuo
       elseif j>ry
           parcial(j) = Theta_EMQ(j)*u(k-j+ry); %Parcial dependente da entrada  
       else
           parcial(j) = Theta_EMQ(j)*y_hat(k-j); %Parcial dependente da saída
       end
       y_hat(k) = y_hat(k)+parcial(j);
   end   
end

MSE = 0;
for i=1:n
    MSE = MSE + i*((y(i)-y_hat(i))^2)/n; %Mean Squared Error
end 

%% Plot
%mostrarEstimativa = 0;
%y_hat = y_hat(1,kmin:n); %Ajustar tamanho para igualar saída real e resíduo

ksi2 = y - y_hat';
var_ksi2 = var(ksi2);
n = length(Y);
if mostrarEstimativa
    t = (0:n-1)';
    hold on;
    plot(t,Y,'DisplayName','Saída Real');
    plot(t,y_hat,'--','DisplayName',['Estimativa ARMAX ' num2str(ry) 'y/',num2str(ru) 'u/', num2str(rn) 'resíduo: ', num2str(MSE) ' MSE']);
    stem(t,ksi2,'Marker','none','DisplayName',['Resíduo da estimativa: variancia=' num2str(var_ksi2)],'LineWidth',1);
    stem(t,ksi,'Marker','none','DisplayName',['Resíduo projetado: iterações=' num2str(iteracoes) ' / variancia=',num2str(var_ksi)],'LineWidth',1);
    legend('off');
    legend('show'); %Update da legenda pra adicionar texto em vez de substituir
    grid on;
end