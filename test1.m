S= [1/9, 1/8, 1/7, 1/6, 1/5, 1/4, 1/3, 1/2, 1, 2, 3, 4, 5, 6, 7, 8, 9];
n = 4;
function P=pairwise_matrix(S,n)

 P=zeros(n,n);

 for i=1:n
 for j=1:n
 if (i==j)
 P(i,i)=1; % fill in the diag
 elseif j>i % fill in the upper triangular array
 P(i,j)=S(ceil(length(S).*rand(1,1)));
 end %eof if

 end % eof for j
 end %eof for i

 %fill in the lower triangular array
 for k=2:n
 for m=1:(k-1)
 P(k,m)=1/P(m,k); %otherwise P(k,m)=inv(P(m,k));
 end %eof k
 end %eof m

endfunction


function [W,CR]=eigenmethod(P)
 [V,lamda]=eig(P);
 [nx,ny]=size(lamda);
 L=zeros(nx,1);

 for i=1:nx
 if imag(lamda(i,i))==0
 L(i)=lamda(i,i);
 end
 end

 [Lmax,ind]=max(L);

 disp(W=abs(V(:,ind)/sum(V(:,ind))));
 n=nx;
 CI=(Lmax-n)/(n-1);

 N=[3 4 5 6 7 8 9 10];
 RI=[0.52 0.89 1.11 1.25 1.25 1.4 1.45 1.49];

 for i=1:length(N)
 if n==N(i)
 CR=CI/RI(i);
 elseif n==1 || n==2
 CR=0;
 end
end %eof for

endfunction



M=15; %number of experts
N=3; %number of criteria
Nalter=3; %number of alternatives
Nf=[2 3 4];% number of factors per criterion
Nfactors=sum(Nf,2);
S=[1/9,1/8,1/7,1/6,1/5,1/4,1/3,1/2,1,2,3,4,5,6,7,8,9]; %nine level scale
Pc=ones(N,N,M); %initialize the PWC matrix of criteria
Pfc1=ones(Nf(1),Nf(1),M); %initialize the PWC matrix of factors under C1
Pfc2=ones(Nf(2),Nf(2),M);%initialize the PWC matrix of factors under C2
Pfc3=ones(Nf(3),Nf(3),M);%initialize the PWC matrix of factors under C3
PAf=ones(Nalter,Nalter,Nfactors,M); %initialize PWC of alternatives relative importance for each factor
%PWC matrix generation
for m=1:M
 Pc(:,:,m)=pairwise_matrix(S,N);
 Pfc1(:,:,m)=pairwise_matrix(S,Nf(1));
 Pfc2(:,:,m)=pairwise_matrix(S,Nf(2));
 Pfc3(:,:,m)=pairwise_matrix(S,Nf(3));

 for i=1:Nfactors
 PAf(:,:,i,m)=pairwise_matrix(S,Nalter);
 end %end of i

end %end of m

w=zeros(N,M);% initialize the matrix of criteria weights - wk (m)
wf1=zeros(Nf(1),M); %initialize the matrix of factor weights under C1 - fj1 (m)
wf2=zeros(Nf(2),M); %initialize the matrix of factor weights under C2 - fj2 (m)
wf3=zeros(Nf(3),M); %initialize the matrix of factor weights under C3 - fj3 (m)
RAf=zeros(Nalter,Nfactors,M); %initialize matrix of alternatives relative importance per factor-Rijk (m)
%Evaluation of weights and relative scores from PWCs
for m=1:M
 w(:,m)=eigenmethod(Pc(:,:,m));
 wf1(:,m)=eigenmethod(Pfc1(:,:,m));
 wf2(:,m)=eigenmethod(Pfc2(:,:,m));
 wf3(:,m)=eigenmethod(Pfc3(:,:,m));

 for i=1:Nfactors
 RAf(:,i,m)=eigenmethod(PAf(:,:,i,m));
 end %end of i

end %end of m
%Estimate average values of weights and relative scores for the M experts
W=mean(w,2); %mean weights of criteria - wk
F1=mean(wf1,2);%mean weights of factors under C1 - fj1
F2=mean(wf2,2);%mean weights of factors under C2 - fj2
F3=mean(wf3,2);%mean weights of factors under C3 - fj3
R=mean(RAf,3); %mean relative scores of alternatives for each factor - Rijk
F=[F1;F2;F3]; %the weights of all the factors in a unique matrix_type
ScenarioValue=zeros(Nalter,1);
for i=1:Nalter %alternatives
 j=0;
 Nfcur=0;

 for k=1:N %criteria
 Nfcur=Nfcur+Nf(k); %Nfcur shows the maximum j for each k

 for j=j+1:Nfcur %Factors
 ScenarioValue(i)=ScenarioValue(i)+W(k)*F(j)*R(i,j);
 end %eof j
 end %eof k

end %eof i

disp('w');
disp(w);
disp('W');
disp(W);
disp('F');
disp(F);
disp('S');
disp(ScenarioValue);
disp('paf');
%disp(PAf);

