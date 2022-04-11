
t1=clock;

[m,n]=size(A);
Q=zeros(m,n);
R=zeros(n,n);
Q2=zeros(m,n);
R2=zeros(n,n);
Q3=zeros(m,n);
R3=zeros(n,n);

%CGS
for k=1:n
    R(1:k-1,k)=Q(:,1:k-1)'*A(:,k);
    v=A(:,k)-Q(:,1:k-1)*R(1:k-1,k);
    R(k,k)=norm(v);
    Q(:,k)=v/R(k,k);
end

%MGS
A2=A;
for k=1:n
    v=A2(:,k);
    R2(k,k)=norm(v);
    Q2(:,k)=v/R2(k,k);
    proj=Q2(:,k)'*A2(:,k+1:n);
    A2(:,k+1:n)=A2(:,k+1:n)-Q2(:,k)*proj;
end

%HS
A3=A;
R3=A3;
H=eye(m);
for k=1:n
    y=R3(k:m,k);
    e=[1;zeros(m-k,1)];
    w=y+y(1,1)/abs(y(1,1))*norm(y)*e;
    v=w/norm(w);
    I=eye(m);
    I(k:m,k:m)=I(k:m,k:m)-2*v*v';
    h=I;
    R3=h*R3;
    H=h*H;
end
Q3=H';

% matlab QR
[Q4,R4]=qr(A);

% 精度
E_CGS=zeros(n,1);
E_MGS=zeros(n,1);
E_HS=zeros(n,1);
E_Matlab=zeros(n,1);

for k=2:n
    E_CGS(k,1)=max(abs(Q(:,1:k-1)'*Q(:,k)));
    E_MGS(k,1)=max(abs(Q2(:,1:k-1)'*Q2(:,k)));
    E_HS(k,1)=max(abs(Q3(:,1:k-1)'*Q3(:,k)));
    E_Matlab(k,1)=max(abs(Q4(:,1:k-1)'*Q4(:,k)));
end

%画图
x=1:1:n;
figure
% 4个方法
% plot(x,E_CGS,'-*r',x,E_MGS,'-*b',x,E_HS,'-*m',x,E_Matlab,'-*k')
% legend('CGS', 'MGS','HS', 'Mat');
% 3个方法
% plot(x,E_MGS,'-*b',x,E_HS,'-*m',x,E_Matlab,'-*k')
% legend('MGS', 'HS', 'Mat');
% %2个方法
plot(x,E_HS,'-*m',x,E_MGS,'-*k')
legend('HS', 'MGS');
% 1个方法
% plot(x,E_HS,'-*b')
% legend('HS');

% det(A)

% %QR求逆：CGS
% CGS_InvB=[];
% for i = 1:size(A,1)
%    b = zeros(1,size(A,1))';
%    b(i)=1;
%    CGS_InvB = [CGS_InvB,backsub(R,Q'*b);];
% end

% %QR求逆：MGS
% MGS_InvB=[];
% for i = 1:size(A,1)
%    b = zeros(1,size(A,1))';
%    b(i)=1;
%    MGS_InvB = [MGS_InvB,backsub(R2,Q2'*b);];
% end

% %QR求逆：HS
% HS_InvB=[];
% for i = 1:size(A,1)
%    b = zeros(1,size(A,1))';
%    b(i)=1;
%    HS_InvB = [HS_InvB,backsub(R3,Q3'*b);];
% end

% %Matlab求逆
% Mat_InvB = inv(A'*A)*A';

% t2=clock;
% etime(t2,t1)

