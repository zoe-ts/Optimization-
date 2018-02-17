function [x_star,f_star] = pso(f,n,max_it,dim)
%An implementation of Particle Swarm Optimization Method  
%
%Inputs: 
%   f: objective function
%   n: number of particles
%   dim: number of dimensions 
%Outputs:
%   x_star: optimal point
%   f_star: optimal value of objective function
%
%Author: Zoe Ts. 2017 

%Initialization
it=1;				%iteration number
y = sym('y', [dim, 1]);
u=zeros(dim,n);
for i=1:n
    for jj=1:dim
        u(jj,i)=0;	%initial velocity
    end
end
d=0.0001;			%deviation
c1=2;				%intelligence constant1
c2=2;				%intelligence constant2
theta_max=0.9;
theta_min=0.4;
%initial positions
lb = -10;
ub = 10;
%initial best value of particles
for i=1:n
	Pbest(1,i)=Inf;	
end
k=0;
%random initial positions of particles
x=zeros(dim,n);
for j=1:n
    for jj=1:dim
        num=(ub-lb).*rand+lb;
        x(jj,j)=num;
    end
end

%PSÏ algorithm
while k<1 && it<max_it

for j=1:n
    %evaluation of f for each particle
	f_eval(it,j)=subs(f,y,x(:,j));
    %best value of each particle
	if f_eval(it,j)<Pbest(1,j)
		Pbest(1,j)=f_eval(it,j);		
		xp_best(:,j)=x(:,j);
	end
end
%global best value of all particles 
[Gbest,Idxbest]=min(Pbest(:));		
xg_best=x(:,Idxbest);
%inertia weight
theta=theta_max-((theta_max-theta_min)/max_it)*it;
r1=rand;
r2=rand;
it=it+1;
%new velocity and position of particles
for j=1:n
	u(:,j)=theta.*u(:,j)+c1*r1.*(xp_best(:,j)-x(:,j))+c2*r2.*(xg_best-x(:,j));		
	x(:,j)=x(:,j)+u(:,j);													
end
%checking of position similarity
for i=1:(n-1)
    if all(abs(x(1,i)-x(1,i+1)))<=d
        k=1;
    end
end

end

x_star=mean(x,dim);
f_star=subs(f,y,x_star);
end