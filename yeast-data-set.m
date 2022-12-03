data=load('yeast.mat');

n=2417;
p=10;

%center the data
datac=data-repmat(sum(data)/n, n, 1);
%find covariant matrix 
covm=cov(datac);

%finding eigenvectors and eigenvalues 
[eigvec, eigval]=eig(covm);

%How many dimensions we should keep? 

%1st method: Scree plot 
%extract diagonal elements 
eigval=diag(eigval);
%order in decending order 
eigval=flipud(eigval);
eigvec=eigvec(:, p:-1:1);

%do a scree plot 
figure, plot(1: length(eigval), eigval, 'ko-')
title('Scree Plot')
xlabel('Eigenvalue Index -k')
ylabel('Eigenvalue')

%2nd method: percentage of variance 
pervar = 100*cumsum(eigval)/sum(eigval);

%3rd broken stick test 
% First get the expected sizes of the eigenvalues.
g = zeros(1,p);
for k = 1:p
for i = k:p
g(k) = g(k) + 1/i;
end
end
g = g/p;
%second step proporsion of variance 
propvar = eigval/sum(eigval);

% Now for the size of the variance.
avgeig = mean(eigval);
% Find the length of ind:
ind = find(eigval > avgeig);
length(ind)

% Using d = 3, we will reduce the dimensionality.
P = eigvec(:,1:3);
Xp = datac*P;
figure,plot3(Xp(:,1),Xp(:,2),Xp(:,3),'k*');
xlabel('PC 1'),ylabel('PC 2'),zlabel('PC 3');



