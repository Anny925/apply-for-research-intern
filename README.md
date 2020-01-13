# apply-for-research-intern
This is a function I wrote for homework in last semester. The function is used to pricing American vanilla options with explicit difference scheme for transformed Black-Scholes PDE. The output c and p are the prices of call and put options respectively, c_delta and p_dejlat are the option delta's for call and put respectively. The input S0 is the initial stock price, X is the strike price, r is the interest rate, T is the terminal time, sigma is volatility, N is number of time step, dx is the step size in x, and [xmin,xmax] gives the truncated domain for the variable x, where x is log of stock price.

function [c,c_delt,p,p_delt]=eds_am_vannilla(S0, X, r, T, sig, N, dx,xmin, xmax)

dt=T/N;
I=(xmax-xmin)/dx;

VGridp=zeros(I+1,N+1);% finite difference grid
VGridc=zeros(I+1,N+1);
ishift=1; 

% Boundary conditions
VGridp(1,:)=max(X-exp(xmin),0); % at S=0;
%put option upper boundary
VGridp(I+1,:)=max(X-exp(xmax),0); % at Smax, unlikely to have positive payoff for large S
VGridc(1,:)=max(exp(xmin)-X,0);
VGridc(I+1,:)=max(exp(xmax)-X,0);

% Terminal condition
VGridp(:,N+1)=max(X-exp((0:I)*dx+xmin),0);
VGridc(:,N+1)=max(exp((0:I)*dx+xmin)-X,0);

alpha=sig^2*dt/(dx^2); smalla=0.5+(r-(sig^2)/2)*dx/2/(sig^2);
c=alpha*smalla/(1+r*dt);
b=(1-alpha)/(1+r*dt);
a=alpha*(1-smalla)/(1+r*dt);

i=(1:I-1)'+ishift;
for n=N:-1:1  % backward time recursive
    VGridp(i,n)=max(a.*VGridp(i-1,n+1)+b.*VGridp(i,n+1)+c.*VGridp(i+1,n+1),X-exp(xmin+(i-ishift)*dx));
    VGridc(i,n)=max(a.*VGridc(i-1,n+1)+b.*VGridc(i,n+1)+c.*VGridc(i+1,n+1),exp(xmin+(i-ishift)*dx)-X);
%    plot(VGrid(i,n)); title(num2str(n));
%    pause(0.15);
end;      
a=(log(S0)-xmin)/dx; a0=floor((log(S0)-xmin)/dx); a1=a0+1;
p0=VGridp(a0+ishift,1); p1=VGridp(a1+ishift,1);
c0=VGridc(a0+ishift,1); c1=VGridc(a1+ishift,1);
p=p0*(a-a1)/(a0-a1)+p1*(a-a0)/(a1-a0);
c=c0*(a-a1)/(a0-a1)+c1*(a-a0)/(a1-a0);
p_delt=(p1-p0)/dx/S0;
c_delt=(c1-c0)/dx/S0;
