% initial matrix
X=[1 0 1 0 0 1 0 0 1 0 0;1 0 0 0 1 1 0 0 0 0 1;0 1 0 1 0 0 1 1 0 0 0;
        1 0 0 0 1 1 0 0 0 1 0;0 1 0 1 0 0 1 0 0 1 0;0 1 0 1 0 1 0 0 1 0 0;
        1 0 1 0 0 0 1 0 0 0 1;1 0 0 0 1 0 1 1 0 0 0;1 0 0 1 0 1 0 1 0 0 0;
        0 1 1 0 0 0 1 0 0 1 0;0 1 0 0 1 1 0 0 1 0 0;1 0 1 0 0 0 1 1 0 0 0];

g =[2 3 2 4] ; % vector of categorical variables
p = length(g); %number of categorical variables
b1 = size(X,2); %number of categories
b2 = sum(g); %input number of categories
if b1 ~= b2;
    disp('Warning:the number of categories in the input vector doesn''t correspond to the given in the indicator matrix. Please, check.');
    return
end
n = size(X,1);
N = sum(sum(X));
P = X/N; %correspondence matrix
r = sum(P,2); %row total (row mass)
c = sum(P,1); %column total (column mass)
Dr = diag(r); %row diagonal matrix
Dc = diag(c); %column diagonal matrix
Z = sqrt(inv(Dr))*(P-r*c)*sqrt(inv(Dc));
%non-trivial eigenvalues (maximum number of dimensions)
if n >= p,
    q = sum(g)-p;
else
    q = p-n;
end
[U,S,V] = svd(Z); %singular values
S;
SV = diag(S); %singular values
SV = SV(1:q); %singular value equal to zero is deleted
EV = SV.^2; %eigenvalues (inertia)
PI = EV./sum(EV)*100; %percent of inertia
PA = cumsum(PI); %cumulative percent of inertia
x2 = N*EV; %chi-square
A = sqrt(Dr)*U;
B = sqrt(Dc)*V;
X2 = N*trace(inv(Dr)*(P-r*c)*inv(Dc)*(P-r*c)'); %total chi-square
v = prod(g-1); %degrees of freedom
% P = 1-chi2cdf(X2,v); %P-value
%disp(' ')
%disp('Inertia and chi-square decomposition of the Multiple Correspondence Analysis')
%disp('for the indicator matrix given.')
%fprintf('------------------------------------------------------------------------\n');
%disp('                                                            Cummulative       ');
%disp(' Singular Value     Inertia      Chi-square     Percent       Percent   ');
%fprintf('------------------------------------------------------------------------\n');
%fprintf(' %10.4f      %10.4f    %10.4f    %10.2f    %10.2f\n',[SV,EV,x2,PI,PA].');
%fprintf('------------------------------------------------------------------------\n');
%fprintf('     Total   %14.4f %13.4f %13.2f\n', sum(EV),N*sum(EV),sum(PI));
%disp(['Variable categories = ',num2str(g)]);
%disp(['                         P-value = ',num2str(P) '; ' 'degrees of freedom = ',num2str(v) '']);
%fprintf('------------------------------------------------------------------------\n');
%disp(' ');
Dm = [];
for i = 1:q
    for s = i+1:q
        if s ~= i
            d = [i s];
            Dm = [Dm;d];
        end
    end
end
co = (q*(q - 1))/2;
gfl = input('Are you interested to get the dimensions plots? (y/n): ','s');

if gfl == 'y'
    disp('--------------');
    fprintf('The pair-wise plots you can get are: %.i\n', co);
    disp('--------------');
    Dm
    fprintf('--------------\n');

    d = 'y';
    while d == 'y';
        p = input('Give me the interested dimensions to plot. Please, use [a b]:');
        figure;
        X = inv(Dr)*A*S; %column coordinates (observations)
        X = X(:,1:q);
        Y = inv(Dc)*B*S'; %row coordinates (categorical variables)
        Y = Y(:,1:q);
        pol = 0;
        hold on
        lg = [];
        for c = 1:length(g),
            plot(Y(1+pol:g(c)+pol,p(1)),Y(1+pol:g(c)+pol,p(2)),'*','Color',paletc(c));
            lg = [lg,['''Levels of categorical variable ' num2str(c) ''',']];
            A1 = Y(1+pol:g(c)+pol,p(1));
            A2 = Y(1+pol:g(c)+pol,p(2));
            for Kl = 1:length(A1),
                text(A1(Kl)+.02,A2(Kl)+.02,num2str(Kl));
            end
            pol = g(c)+pol;
        end
        lg(end) = ' ';
        eval(['legend(' lg ')']);
        title('Plot of variable levels of the performed Multiple Correspondence Analysis.');
        xlabel(['Dimension ' num2str(p(1))';]);
        ylabel(['Dimension ' num2str(p(2))';]);

        disp(' ')
        vl = input('Are you interested to plot the origin reference of the variable levels? (y/n): ','s');
        if vl == 'y'
            disp(' ')
            disp('Note: At the end of the program execution. If you are interested to fit the variable labels')
            disp('on the generated figures you can turn-on active button ''Edit Plot'', do click on the')
            disp('selected label and drag to fix it on the desired position. Then turn-off active ''Edit Plot''.')
            disp(' ')
            pol = 0;
            for c = 1:length(g),
                hold on
                A1 = Y(1+pol:g(c)+pol,p(1));
                A2 = Y(1+pol:g(c)+pol,p(2));
                plot2org(A1,A2,'--k'), %this m-file was taken from Jos's plot2org
                pol = g(c)+pol;
            end
        else
        end
        o = input('Are you interested to plot the origin quadrature? (y/n): ','s');
        if o == 'y'
            hline(0,'k'); %this m-file was taken from Brandon Kuczenski's hline and vline zip
            vline(0,'k'); %this m-file was taken from Brandon Kuczenski's hline and vline zip
            disp(' ')
            d=input('Do you need another plot? (y/n):','s');
        else
        end
    end
else
end
function res = paletc(n)
% - utility program for selecting plot colors
while (n>10)
   n = n - 10;
end
% build cell-array of colors to choose from:
colors{1}  = [0 0 1];          % blue
colors{2}  = [1 0 0];          % red
colors{3}  = [0 1 0];          % green
colors{4}  = [0.5 0 0.5];      % purple
colors{5}  = [0 0.5 0];        % dark green
colors{6}  = [0.28 0.73 0.94]; % sky
colors{7}  = [0 0 0];          % black
colors{8}  = [0 0 0.5];        % navy
colors{9}  = [1 0.5 0];        % orange
colors{10} = [1 0 1];          % magenta
colors{11} = [1 1 0];          % yellow
res = colors{n};
end
