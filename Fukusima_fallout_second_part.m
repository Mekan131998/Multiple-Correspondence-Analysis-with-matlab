format long
r=log(2)/(8); !half-life of I-131 is 8 days 

dt=5/100;
w=1/5;
f=1/(3*365);
h=1/10;
g=1/15;
p=1/7;
u=1/30;
x=1/(3*30);
y=1/(4*30);
z=1/16;
e=1/2;

n=4000; %steps to be plot, M^22000

Hf=zeros(1,n);
Hp=zeros(1,n);


M=[1 0 0 0 0 0 0 0 0 z*dt 0 0 
   0 1 0 0 0 0 0 x*dt 0 0 0 0 
   0 0 1 0 0 0 0 0 0 0 0 e*dt 
   0 0 0 1 0 0 0 0 0 0 0 r*dt
   0 0 0 0 1 r*dt r*dt r*dt r*dt r*dt r*dt 0
   0 0 0 0 0 1-(y+g+w+r)*dt 0 0 0 0 0 0
   0 0 0 0 0 w*dt 1-(f+r)*dt 0 0 0 0 0 
   0 0 0 0 0 0 f*dt 1-(x+r)*dt 0 0 0 0 
   0 0 0 0 0 g*dt 0 0 1-(p+u+r)*dt 0 0 0
   0 0 0 0 0 0 0 0 p*dt 1-(z+r+u)*dt 0 0
   0 0 0 0 0 0 0 0 u*dt u*dt 1-r*dt 0 
   0 0 0 0 0 y*dt 0 0 0 0 0 1-(e+r)*dt];
   

N=zeros(12,n+1);
N(6,1)=1;
for j=2:n+1
    N(:,j)=M*N(:, j-1);
    Hf(1,j)=Hf(1,j-1)+x*dt*N(6,j-1);
    Hp(1,j)=Hp(1,j-1)+y*dt*N(7,j-1);
end

figure(1)
for i=1:12
    plot(N(i,:));
    legend("Stable", "Human", "Air","Cow", "Ground", "Vegetable");
    hold on;
end

x_axis = 1:n+1;
figure(2)
plot (x_axis, Hf,'r-', x_axis, Hp,'g-')
legend("H_f", "H_p");
