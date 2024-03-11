function hess=hess_svm1(x,y,theta,Mchol_1,Mchol_2)
% Hessian of the objective Function
% Support vector machine with Cobb-Douglas function
% f(w,b,ka_1,ka_2)=-theta*log(ka_1^2/(ka_1^2+1))-(1-theta)*log(ka_2^2/(ka_2^2+1))
% Input: 
%          x=(w,b,ka_1,ka_2) in R^{n+3}
%          theta := parameter of model
C1=theta;
C2=1-theta;
x=x(:);
n=length(x);

z1=2*(3*x(n-1)^2+1)/((x(n-1)^2+1)*x(n-1))^2;
z2=2*(3*x(n)^2+1)/((x(n)^2+1)*x(n))^2;
Hf=sparse(n,n);
Hf(n-2,n-2)=0;
Hf(n-1,n-1)= C1*z1;
Hf(n,n)= C2*z2;


Hg11= - sparse(n,n)*y(1);
Hg12=sparse(n,n);
for j=1:(n-3)
    Hg12(n-2,j) = - Mchol_1(j,j)*y(j+1);
    Hg12(j,n-2) = Hg12(n-2,j);
end
Hg21= - sparse(n,n)*y(n-1);
Hg22=sparse(n,n);
for j=1:(n-3)
    Hg22(n-2,j) = - Mchol_2(j,j)*y(n-1+j);
    Hg22(j,n-2) = Hg22(n-2,j);
end
Hg3=-sparse(n,n)*y(2*n-3);
Hg4=-sparse(n,n)*y(2*n-2);

hess = Hf +Hg11+ Hg12+Hg21+Hg22+Hg3+Hg4;

%modified newton step
epsilon = min(eig(hess));
if epsilon <= 0.0
    hess = hess+(abs(epsilon)+0.1)*eye(n);
end

return
