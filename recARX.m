function [MSE,y_hat,Psi,Theta,ksi] = recARX(y,u,ry,ru,Ts,mostrarEstimativa)

%Numero de Amostras
n = length(u);

%kmin
if ru>=ry
    kmin = ru+1;  
else
    kmin = ry+1;
end

%Numero de Regressores
rtotal = ry+ru;

%Valor inicial dos parametros
Theta = zeros(rtotal,kmin-1);
%Valor inicial da Matriz P
P(:,:,kmin-1) = eye(rtotal);

%Pre alocacoes
Psi = zeros(rtotal,n-kmin);
%Loop principal do Estimador MQR
for k=kmin:n   
     %Vetor de regressao Psi
     %Parcial de Psi dependente da saida
     Psiy = zeros(1,ry); %Pre alocacao
     for i=0:ry-1
        Psiy(1,i+1) = y(k-i,1); 
     end
     %Parcial de Psi dependente da entrada
     Psiu = zeros(1,ru); %Pre alocacao
     for i=0:ru-1
        Psiu(1,i+1) = u(k-i,1); 
     end
     Psi(:,k) = [Psiy Psiu]'; %Vetor de Regressao
     
     %Equacoes:
     % Eq1 -> ganho(k) = P(k-1) . psi(k) . [1 + psi(k)' . P(k-1) . psi(k)]^-1 
     % Eq2 -> inovacao(k) = y(k) - psi(k)' . Theta(k-1)
     % Eq3 -> Theta(k) = Theta(k-1) + ganho(k) . inovacao(k)
     % Eq4 -> P(k) = P(k-1) - ganho(k) . psi(k)' . P(k-1)
   
     %Eq1 Ganho dos MQR
     ganho(:,k) = P(:,:,k-1)*Psi(:,k) * [1 + Psi(:,k)'*P(:,:,k-1)*Psi(:,k)]^(-1);
     %Eq2 Inovacao
     inov(k) = y(k) - Psi(:,k)'*Theta(:,k-1);
     %Eq3 Atualizacao dos Parametros do modelo
     Theta(:,k) = Theta(:,k-1) + ganho(:,k)*inov(k);
     %Eq4 Atualizacao da matriz P (Covariancia)
     P(:,:,k) = P(:,:,k-1) - ganho(:,k)*Psi(:,k)'*P(:,:,k-1);
     
  %SINAL DEVE SER NEGATIVO -
     
     %Estimativa da saida atual
     y_hat(k) = Psi(:,k)'*Theta(:,k);
end

ksi = y' - y_hat;

MSE = 0;
for i=1:n
    MSE = MSE + i*((y(i)-y_hat(i))^2)/n; %Mean Squared Error
end 

if mostrarEstimativa
    %tempo
    t = [0:Ts:(n-1)*Ts]';
    
    figure();
    subplot(311);
    for i=1:rtotal
        plot(t,Theta(i,:)','DisplayName',['Theta ' num2str(i)]);
        legend('off');
        legend('show'); %Update da legenda pra adicionar texto em vez de substituir
        hold on;
    end
    title(['parametros: ' num2str(ry) 'y/' num2str(ru) 'u']);
    ylabel('amplitde');
    xlabel('tempo [s]');

    subplot(312);
    plot(t,inov);
    title('inovacao');
    ylabel('amplitde');
    xlabel('tempo [s]');

    subplot(313);
    plot(t,ganho(1,:),'r');
    title('ganho');
    ylabel('ganho');
    xlabel('tempo [s]');

    figure();
    plot(t,y,'r');
    hold on;
    plot(t,y_hat','-.b');
    legend('Dado',['Estimativa ' num2str(ry) 'y/' num2str(ru) 'u  MSE=' num2str(MSE)]);
    grid on;
end