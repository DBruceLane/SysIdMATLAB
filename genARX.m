function [fitness,MSE,y_hat,Psi,Theta,ksi] = genARX(y,u,ry,ru,mostrarEstimativa)
% Generalizacao de modelagem ARX
n = size(u,1);

% Se o menor indice na equacao de y(k) eh (k-1), o menor k possível eh 2 (pro Matlab iniciar indice em 1)
% Se fosse o caso kmin = 2

if ru>=(ry+2)
    kmin = ru;  
    %Quando ha 2 regressores a mais de entrada que saida, ex: ry = 1 e ru = 3
    %y(k) = y(k-1) + u(k) + u(k-1) + u(k-2) 
    %Nesse caso observa-se o ultimo elemento de entrada, kmin seria 3
else
    kmin = 1+ry;
    %Nos outros casos observa-se o ultimo elemento de saida
end

Y = y(kmin:n);

Psiy = zeros(1+n-kmin,ry); %Elementos de saida na matriz Psi
for i=1:ry
        Psiy(1:1+n-kmin,i) = y(kmin-i:n-i); 
end

Psiu = zeros(1+n-kmin,ru); %Elementos de entrada na matriz Psi
for i=1:ru
       Psiu(1:1+n-kmin,i) = u(kmin-i+1:n-i+1);     
end

Psi = [Psiy Psiu];
Theta = inv(Psi'*Psi)*Psi'*Y; %Pseudo Inversa

%Free Run Simulation
y_hat = zeros(1,ry);
for i=1:ry
    y_hat(i) = y(i); %Condicoes de Contorno
end

parcial = zeros(1,ry+ru); %Parciais de y_hat
%Exemplo de equacao pra mostrar os indices do for
%y_hat(k) = Theta(j=1)*y(k-1) + Theta(j=2)*y(k-2) + Theta(j=3)*u(k) + Theta(j=4)*u(k-1)
for k=kmin:n  
   y_hat(k) = 0;
   for j=1:ry+ru
       if j>ry
           parcial(j) = Theta(j)*u(k-j+ry+1); %Parcial dependente da entrada  
       else
           parcial(j) = Theta(j)*y_hat(k-j); %Parcial dependente da saida
       end
       y_hat(k) = y_hat(k)+parcial(j);
   end   
end

MSE = 0;
for i=1:n
    MSE = MSE + i*((y(i)-y_hat(i))^2)/n; %Mean Squared Error
end 

ksi = y - y_hat';

NRMSE = goodnessOfFit(y_hat',y,'NRMSE');
fitness = 100*(1-NRMSE);

if mostrarEstimativa
    t = (0:n-1)';
    hold on;
    plot(t,y_hat,'--','DisplayName',[num2str(ry) 'y/',num2str(ru) 'u: ',num2str(MSE) ' MSE']);
    stem(ksi,'Marker','none','DisplayName','Residuo');
    legend('off');
    legend('show'); %Update da legenda pra adicionar texto em vez de substituir
    grid on;
end
