format long
r=log(2)/8; ! a half-life of iodine-131 constitutes 8 days
w=1/5;
f=1/(3*365); !1/year
h=1/10;
g=1/15;
p=1/7;
u=1/(1*30);!1/month
x=1/(3*30);!1/month
y=1/(4*30);!1/month
z=1/16;
e=1/2;

dt=5/100;

n=4000; !steps to be plot, M^22000

Hp=zeros(1,n);
Hf=zeros(1,n);

M=[1 0 0 0 0 0 0 0 0 e*dt 
   0 1 0 0 0 0 0 0 0 r*dt
   0 0 1 r*dt r*dt r*dt r*dt r*dt r*dt 0
   0 0 0 1-(y+g+w+r)*dt 0 0 0 0 0 0
   0 0 0 w*dt 1-(f+r)*dt 0 0 0 0 0
   0 0 0 0 f*dt 1-(x+r)*dt 0 0 0 0
   0 0 0 g*dt 0 0 1-(p+u+r)*dt 0 0 0
   0 0 0 0 0 0 p*dt 1-(r+u)*dt 0 0
   0 0 0 0 0 0 u*dt u*dt 1-r*dt 0
   0 0 0 y*dt 0 x*dt 0 0 0 1-(e+r)*dt];

N=zeros(10,n+1);
N(5,1)=1;
for j=2:n+1
    N(:,j)=M*N(:, j-1);
    Hp(1,j)=Hp(1,j-1)+x*dt*N(5,j-1);
    Hf(1,j)=Hf(1,j-1)+y*dt*N(6,j-1);
end

figure(1)
for i=1:9
    plot(N(i,:));
    legend("Sea Water", "Fish", "Ground", "Plants", "Underground", "Human", "Air", "Urine", "Stable");
    hold on;
end

x_axis = 1:n+1;
figure(2)
plot (x_axis, Hp,'g-', x_axis, Hf,'r-')
legend("H_p", "H_f");
