function W = OLSR(fea,Y)
% fea 每一列是一个样本
% gnd 类标签，随后转换为类指示矩阵
% d 样本最后的维数
% Y = indicatedmatrix(gnd);
[n,m] = size(fea); %m是样本的维数
[~,c] = size(Y);
X = fea';

%--------------------------%
H = eye(n,n)-(1/n)*ones(n,n);


%------EB算法求W-----------%
[V,D] = eig((H*X')'*(H*X')); %最小n-m个特征值对应的特征向量
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

