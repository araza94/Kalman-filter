close all
clear all
load KF_dataset_1
x = zeros(length(t)+1,3);
P = zeros(3,3,length(t)+1);
P(:,:,1) = Po;
y_est = zeros(length(t),1);
Py = zeros(length(t),1);
for k=2:length(t)+1
    x(k,:) = Phi*x(k-1,:)' + Gamma*u(k-1);
    P(:,:,k) = Phi*P(:,:,k-1)*Phi' + Q*eye(3);
    Py(k-1) = H*P(:,:,k)*H' + R;
    K = P(:,:,k)*H'/(H*P(:,:,k)*H'+R);
    x(k,:) = x(k,:)' + K*(y(k-1) - H*x(k,:)');
    P(:,:,k) = P(:,:,k) - K*H*P(:,:,k);
    y_est(k-1) = H*x(k,:)';
end
r=y-y_est;
z=r./sqrt(Py);

figure;
plot(t,x(2:402,1),'k-',t,x(2:402,1)+sqrt(reshape(P(1,1,2:402),[],1)),'k--',t,x(2:402,1)-sqrt(reshape(P(1,1,2:402),[],1)),'k:')
xlabel('t')
legend('x_{1}','x_{1}+\sigma_{x_{1}}','x_{1}-\sigma_{x_{1}}')
figure;
plot(t,x(2:402,2),'k-',t,x(2:402,2)+sqrt(reshape(P(2,2,2:402),[],1)),'k--',t,x(2:402,2)-sqrt(reshape(P(2,2,2:402),[],1)),'k:')
xlabel('t')
legend('x_{2}','x_{2}+\sigma_{x_{2}}','x_{2}-\sigma_{x_{2}}')
figure;
plot(t,x(2:402,3),'k-',t,x(2:402,3)+sqrt(reshape(P(3,3,2:402),[],1)),'k--',t,x(2:402,3)-sqrt(reshape(P(3,3,2:402),[],1)),'k:')
xlabel('t')
legend('x_{3}','x_{3}+\sigma_{x_{3}}','x_{3}-\sigma_{x_{3}}')

figure;
plot(t,y,'k-.',t,y_est,'k-',t,y_est+sqrt(Py),'k--',t,y_est-sqrt(Py),'k:')
xlabel('t')
legend('y','')
legend({'$y$','$\hat{y}$','$\hat{y}+\sigma_{y}$','$\hat{y}-\sigma_{y}$'})
set(legend,'Interpreter','latex')

figure;
plot(t,r,'k-',t,r+sqrt(Py),'k--',t,r-sqrt(Py),'k:')
xlabel('t')
legend('r','r+\sigma_{r}','r-\sigma_{r}')

figure;
plot(t,z,'k-')
xlabel('t')
legend('r/\sigma_{r}')
