%% This is the setup code. It takes a starting vertex and generates N2, N3, and that matrix
coldict=["a1","c1","e1","g1","b2","d2","f2","h2","a3","c3","e3","g3","b4","d4","f4","h4","a5","c5","e5","g5","b6","d6","f6","h6","a7","c7","e7","g7","b8","d8","f8","h8"
];
rowdict=["b1","d1","f1","h1","a2","c2","e2","g2","b3","d3","f3","h3","a4","c4","e4","g4","b5","d5","f5","h5","a6","c6","e6","g6","b7","d7","f7","h7","a8","c8","e8","g8"
];
% This finds N1 for a white starting vertex (in a row) Doing this for a black starting vertex would just require us to change the rows and columns, very easy
N1=zeros(1,0);
A=table2array(KADJ);
N0=1; %1 is b1 here
for c = 1:32
    if A(N0,c)==1
        N1(end+1)=c;
    end
end
N1;
% This finds N2
N2=zeros(1,0);
for c=N1
    for r = 1:32
        if A(r,c)==1
        N2(end+1)=r;
        end
    end
end
N2=setdiff(N2,N0);
% This finds N3
N3=zeros(1,0);
for r=N2
    for c = 1:32
        if A(r,c)==1
        N3(end+1)=c;
        end
    end
end
N3=setdiff(N3,N1);
% This creates the matrix for the N2 N3 graph by making rows of zeros for the other vertices. 
% This in effect makes a matrix for the graph where the vertices remain,
% but the edges are just removed.
T=1:32;
N2N3=A;
N2N3(setdiff(T,N2),:)=0;
N2N3(:,setdiff(T,N3))=0;       
%% Test for coverage by size n
for n=1:10
v = 1:length(N3);
Dn = nchoosek(v,n);
covn=zeros(length(Dn(:,1)),1);
cov=zeros(length(N2),1);
%This creates the matrix from the N2 N3 graph with only columns at the defenders and rows  N2. 
T=1:32;
for r=1:length(Dn(:,1))
    DEFMAT=N2N3;
    DEFMAT(setdiff(T,N2),:)=[];
    DEFMAT(:,setdiff(T,Dn(r,:)))=[]; 
    cov=zeros(length(N2),1);
        for r=1:length(N2)
        cov(r,1)=sum(DEFMAT(r,:));
        end
    covn(r,1)=all(cov);
end
if nnz(covn)~=0
    break
end
end
GC=find(~all(covn==0,2));
for t=GC
    GCn=[GCn;Dn(t,:)];
end
    for i=1:length(GCn(:,1))
a(i,:)=coldict(GCn(i,:));
    end
a