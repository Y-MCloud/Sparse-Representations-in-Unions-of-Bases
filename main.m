clc;
clear;
n = 4;              %the exponent of q
q = 2^n;            %the order q
GF = gf(0:q-1,n);   %F_q: Galois Field of order q
D = [];             %the dictionary D
x = [];             %the sparse vector x


%The addition table for F_q
add_mat=[];
for i=1:q
    for j=1:q
        tempadd=GF(i)+GF(j);
        add_mat(i,j)=double(tempadd.x);
    end
end

%The multiplication table for F_q
mult_mat=GF'*GF;
dmult_mat=double(mult_mat.x);

%Permutated Hadamard matrix of order q
Hadamard = hadamard(q)/sqrt(q);
Hadamard(2:q,:)=Hadamard(q:-1:2,:);

%construct the dictionary D
for i=1:q
    for j=1:q
        col_LatinSquare=mult_mat(:,i)+GF(j);
        LatinSquare(:,j)=double(col_LatinSquare.x); %Latin Square over F_q
    end
    B = constructBase(LatinSquare+1,Hadamard,q); %the corresponding base of LatinSquare and Hadamard
    D = [D,B];  %unite the bases B_0,B_1,...,B_{q-1}
end
I = eye(q);
B_infinity = kron(I,Hadamard);  %the bases B_{\infinity}
D = [D,B_infinity];   %unite the bases B_{\infinity}

%construct the sparse vector x
index=diag(dmult_mat)+1;
for i=1:q
    x=[x;kron(I(:,index(i)),I(:,i))];
end
x=[x;-kron(I(:,1),I(:,1))];

Dx=D*x;






