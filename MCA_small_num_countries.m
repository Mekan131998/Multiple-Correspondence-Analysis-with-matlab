%Our data matrix for small number of countries 
N=[
1	0	1	0	0	1	0	1	0	1	0	0	1	1	0	0	0	0	1	0	0	1	0	0;
1	0	0	1	0	1	1	0	0	1	0	0	1	1	0	0	0	0	0	0	1	0	0	1;
0	1	0	1	1	0	0	1	1	0	0	1	0	0	0	1	0	0	0	1	0	0	1	0;
1	0	1	0	0	1	0	1	1	0	1	0	0	0	0	1	1	0	0	0	0	1	0	0;
1	0	0	1	0	1	0	1	1	0	0	1	0	0	1	0	0	1	0	0	0	0	1	0;
    ];

%adding columns:
n_i=sum(N,2);

%adding rows:
n_j=sum(N);

[nc, nr]=size(N);
%sum
n=sum(sum(N));

%Uncorrolated matrix
U=(n_i*n_j)/n;

%sizes
[numCols,numRows] = size(N);

%X matrix
X=zeros(numCols,numRows);
for i=1:numCols
    for j=1:numRows
        X(i,j)=X(i,j)+(N(i,j)-U(i,j))/sqrt(U(i,j));
    end
end

%centering matrix:
[n_x, p_x]=size(X);
xbar=mean(X);

Xc=X-repmat(xbar, n_x, 1);

% PCA of Xc is CA of N

%V vector
V=(Xc'*Xc)/(n-1);

%eigenvalues
[E, D]=eig(V);

%diagonal matrix
Dh=diag(n_i);
Df=diag(n_j);
%rows in P1P2
R=sqrt(n)*inv(Dh)*N*sqrt(inv(Df))*E;
C=sqrt(n-1)*sqrt(inv(Df))*E*sqrt(D);

%1st method: Scree plot
%extract diagonal elements
eigval=diag(D);
%order in decending order
eigval=flipud(eigval);
eigvec=E(:, nr:-1:1);

%do a scree plot
figure, plot(1: length(eigval), eigval, 'ko-')
title('Scree Plot')
xlabel('Eigenvalue Index -k')
ylabel('Eigenvalue')

%2nd method
propvar=eigval/sum(eigval);

% 3rd method:
avgeig = mean(eigval);
% Find the length of ind:
ind = find(eigval > avgeig);
length(ind);

P=E(:,1:8);
%P1, P2, P3, P4, P5

Xp=Xc*P;

[n_Xp, p_Xp]=size(Xp);

%graph of dimension 1 vs dimension 2
figure, plot(Xp(:,1), Xp(:,2),'b*')
xlabel('PCA1 39.56%'), ylabel('PCA2 27.80%')
for i=1:n_Xp
        text(Xp(i,1),Xp(i,2),num2str(i));
end

%graph of dimension 2 vs dimention 3
figure, plot(Xp(:,2), Xp(:,3),'g*')
xlabel('PC2 27.80%'), ylabel('PC3 18.89%')
for i=1:n_Xp
        text(Xp(i,2),Xp(i,3),num2str(i));
end

%graph of dimension 2 vs dimention 3
figure, plot(Xp(:,1), Xp(:,3),'g*')
xlabel('PCA1 39.56%'), ylabel('PC3 18.89%')
for i=1:n_Xp
        text(Xp(i,1),Xp(i,3),num2str(i));
end
