function W = OLSR(fea,Y)
% fea ÿһ����һ������
% gnd ���ǩ�����ת��Ϊ��ָʾ����
% d ��������ά��
% Y = indicatedmatrix(gnd);
[n,m] = size(fea); %m��������ά��
[~,c] = size(Y);
X = fea';

%--------------------------%
H = eye(n,n)-(1/n)*ones(n,n);


%------EB�㷨��W-----------%
[V,D] = eig((H*X')'*(H*X')); %��Сn-m������ֵ��Ӧ����������
D = diag(D); 
[ ~,I] = sort(D,'ascend');
E = V(:,I(1:m-c));
iter=1;
flag=1;
while flag==1;
iter=iter+1;
[U,~,V] = svd((H*X')'*[H*Y,(H*X')*E]);
G = U*V';
W = G(:,1:c);
E = G(:,(c+1:m));
obj(1,iter) = norm(H*X'*W-H*Y,'fro')^2;
if (abs(obj(1,iter)-obj(1,iter-1))<10e-5)
flag=0;
end
end

