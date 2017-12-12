function [x_star,f_star] = univariate(f,n,x,x_k)
%An implementation of Univariate Method for optimization 
%
%Inputs: 
%   f: objective function
%   n: number of dimensions
%   x: symbolic vector of n size
%   x_k: initial point
%Outputs:flet
%   x_star: optimal point
%   f_star: optimal value of objective function
%
%Author: Zoe Ts. 2017 

% Initialization
e=0.01;			%minimum allowable step length 
i=1;			%iteration
s=[1;0];		%initial search direction

%calculating starting f values
f_new=subs(f,x,x_k);				%calculating f
x_kplus=x_k+e*s;
fplus=subs(f,x,x_kplus);			%calculating f+						
x_kminus=x_k-e*s;	
fminus= subs(f,x,x_kminus);			%calculating f-

%correct direction s and step 
syms l;
if fplus<f_new
	g= subs(f, x, (x_k+s*l));
	newg=diff(g,l)==0;
	optl=solve(newg,l);		%optimal step length
	x_k1=x_k+optl*s;
	f_new=subs(f,x,x_k1);
elseif fminus<f_new
	g= subs(f, x, (x_k-s*l));
	newg=diff(g,l)==0;
	optl=solve(newg,l);		%optimal step length
	x_k1=x_k-optl*s;
	f_new=subs(f,x,x_k1);
end
x_k=x_k1;
i=i+1;

f_pr=10;	
jj=1;
% Univariate algorithm
while f_pr>f_new
	%calculating search direction
	if i<=n
		j=i;
		s=zeros(n,1);
		s(j,1)=1;					
	elseif i>n
		s=zeros(n,1);
		s(jj,1)=1;
		jj=jj+1;
		if jj>n
			jj=1;
		end
	end
	f_new=subs(f,x,x_k);				%calculating f
	f_pr=f_new;
	x_kplus=x_k+e*s;
	fplus=subs(f,x,x_kplus);			%calculating f+	
	x_kminus=x_k-e*s;
	fminus=subs(f,x,x_kminus);			%calculating f-	
	%correct direction s and step 
	syms l;
    if fplus<f_new
		g= subs(f, x, (x_k+s*l));
		newg=diff(g,l)==0;
		optl=solve(newg,l);		%optimal step length
		x_k1=x_k+optl*s;
		f_new=subs(f,x,x_k1);	%calculating new f
	elseif fminus<f_new
		g= subs(f, x,(x_k-s*l));
		newg=diff(g,l)==0;
		optl=solve(newg,l);		%optimal step length
		x_k1=x_k-optl*s;
		f_new=subs(f,x,x_k1);	%calculating new f
    end
	x_k=x_k1;
	i=i+1;
end

x_star=x_k;
f_star=f_new;
end







