%%prior to running this code in Matlab the file KADJ.xlsx must be loaded

coldict=["a1","c1","e1","g1","b2","d2","f2","h2","a3","c3","e3","g3","b4","d4","f4","h4","a5","c5","e5","g5","b6","d6","f6","h6","a7","c7","e7","g7","b8","d8","f8","h8"
];
rowdict=["b1","d1","f1","h1","a2","c2","e2","g2","b3","d3","f3","h3","a4","c4","e4","g4","b5","d5","f5","h5","a6","c6","e6","g6","b7","d7","f7","h7","a8","c8","e8","g8"
];

%% clears all variables so you can rerun
clear N1 c n1a n2a n3a N2N3 N3 p r s def ALG s N2 N0 k mindegvert degN2 vdef vdefn T remove
%% This finds N1 for a white starting vertex (in a row) Doing this for a black starting vertex would just require us to change the rows and columns, very easy
N1=zeros(1,0);
A=transpose(KADJ);
N0=13; %1 is b1 here
for c = 1:32
    if A(N0,c)==1
        N1(end+1)=c;
    end
end
for i=1:length(N1(:,1))
n1a(i,:)=coldict(N1(i,:));
end
n1a
%% This finds N2
N2=zeros(1,0);
for c=N1
    for r = 1:32
        if A(r,c)==1
        N2(end+1)=r;
        end
    end
end
N2=setdiff(N2,N0);
for i=1:length(N2(:,1))
n2a(i,:)=rowdict(N2(i,:));
end
n2a
%% This finds N3
N3=zeros(1,0);
for r=N2
    for c = 1:32
        if A(r,c)==1
        N3(end+1)=c;
        end
    end
end
N3=setdiff(N3,N1);
for i=1:length(N3(:,1))
n3a(i,:)=coldict(N3(i,:));
end
n3a
%% This creates the matrix for the N2 N3 graph by making rows of zeros for the other vertices. 
% This in effect makes a matrix for the graph where the vertices remain,
% but the edges are just removed.
T=1:32;
N2N3=A;
N2N3(setdiff(T,N2),:)=0;
N2N3(:,setdiff(T,N3))=0;        

%% N2 is rows N3 is columns. 
ALG=N2N3;
def=zeros(1,0);

%% Compute the degree of vertices in N2 (rows) 
degN2=zeros(32,1);
mindegvert=zeros(1,0);
for i=1:32
    degN2(i)=sum(ALG(i,:));
end
%creates a vector of the N2 vertices with minimal degree
for i=1:32
    if degN2(i)==min(degN2(degN2>0))
        mindegvert(end+1)=i;
    end
end
    sum(degN2)
%% Step one of Algorithm: Identify all forced defenders. 
for r=1:32
    if sum(ALG(r,:))==1
        %if forced defender, record this column as a defender
        c=find(ALG(r,:)==1);
        def(end+1)=c;
        %replace this column with zeros and replace all rows defended by it
        %with zeros
        %This finds all rows defended by the forced defender and zeros the
        %row out
        for p=find(ALG(:,c)==1)
            ALG(p,:)=0;
        end
    end
end


%% remove N3 of degree 1
remove=zeros(1,0);
for c=1:32
    if sum(ALG(:,c))==1
        ALG(:,c)=0;
        remove(end+1)=c;
    end
end

%% Recheck for forced defenders. 
for r=1:32
    if sum(ALG(r,:))==1
        %if forced defender,co record this column as a defender
        c=find(ALG(r,:)==1);
        def(end+1)=c;
        %replace this column with zeros and replace all rows defended by it
        %with zeros
        %This finds all rows defended by the forced defender and zeros the
        %row out
        for p=find(ALG(:,c)==1)
            ALG(p,:)=0;
        end
    end
end
%% Compute the degree of vertices in N2 (rows) 
degN2=zeros(32,1);
mindegvert=zeros(1,0);
for i=1:32
    degN2(i)=sum(ALG(i,:));
end
%creates a vector of the N2 vertices with minimal degree
for i=1:32
    if degN2(i)==min(degN2(degN2>0))
        mindegvert(end+1)=i;
    end
end
        
%% check for neighbor degree of minimum degree vertices in N2

for k=mindegvert
    vdef=zeros(1,0);
    vdefn=zeros(1,0);
    for i=1:32
        if ALG(k,i)==1
            vdef(1,end+1)=i;
            vdefn(1,end+1)=sum(ALG(:,i));
        end
    end
end



%% choose defender

%this identifies the defenders with maximum vertex degree. If there's a tie
%it creates a vector with all the tied N3 vertices. It then chooses the
%defender by looking at the maximum number of the index of the defender
%(fix this to look at neighbor degrees)
c=vdef(find(vdefn==max(vdefn)));
d=max(c);
def(end+1)=d;

%this deletes everything defended by the new defender by deleting it's row
for p=find(ALG(:,d)==1)
            ALG(p,:)=0;
end

%% Test for coverage

defend=[13,14,18,11];
DEFMAT=N2N3;
%This creates the matrix from the N2 N3 graph with only columns at the defenders and rows  N2. 
T=1:32;
DEFMAT(setdiff(T,N2),:)=[];
DEFMAT(:,setdiff(T,defend))=[];        

cov=zeros(length(N2),1);
for r=1:length(N2)
    cov(r,1)=sum(DEFMAT(r,:));
end
all(cov)

%% Test for coverage by size 3
v = 1:length(N3);
D3 = nchoosek(v,4);
cov3=zeros(length(D3(:,1)),1);
cov=zeros(length(N2),1);
%This creates the matrix from the N2 N3 graph with only columns at the defenders and rows  N2. 
T=1:32;
for r=1:length(D3(:,1))
    DEFMAT=N2N3;
    DEFMAT(setdiff(T,N2),:)=[];
    DEFMAT(:,setdiff(T,D3(r,:)))=[]; 
    cov=zeros(length(N2),1);
        for s=1:length(N2)
        cov(s,1)=sum(DEFMAT(s,:));
        end
    cov3(r,1)=all(cov);
end
any(cov3)
%zero is no coverage, one is coverage
%% Test for coverage by size 3
v = 1:length(N3);
D3 = nchoosek(v,4);
cov3=zeros(length(D3(:,1)),1);
cov=zeros(length(N2),1);
%This creates the matrix from the N2 N3 graph with only columns at the defenders and rows  N2. 
T=1:32;
for r=1:length(D3(:,1))
    DEFMAT=N2N3;
    DEFMAT(setdiff(T,N2),:)=[];
    DEFMAT(:,setdiff(T,D3(r,:)))=[]; 
    cov=zeros(length(N2),1);
        for r=1:length(N2)
        cov(r,1)=sum(DEFMAT(r,:));
        end
    cov3(r,1)=all(cov);
end
nnz(cov3)
%zero is no coverage, one is coverage
%% Test for coverage by size 4
v = 1:length(N3);
D4 = nchoosek(v,4);
cov4=zeros(length(D4(:,1)),1);
cov=zeros(length(N2),1);
GC4=zeros(0,0);
%This creates the matrix from the N2 N3 graph with only columns at the defenders and rows  N2. 
T=1:32;
for r=1:length(D4(:,1))
    DEFMAT=N2N3;
    DEFMAT(setdiff(T,N2),:)=[];
    DEFMAT(:,setdiff(T,D4(r,:)))=[]; 
    cov=zeros(length(N2),1);
        for s=1:length(N2)
        cov(s,1)=sum(DEFMAT(s,:));
        end
    cov4(r,1)=all(cov);
end
nnz(cov4)
GC=find(~all(cov4==0,2));
for t=GC
    GC4=[GC4;D4(t,:)];
end
GC4
    for i=1:length(GC4(:,1))
C4t(i,:)=rowdict(GC4(i,:));
    end
C4t
%% Test for coverage by size 5
v = 1:length(N3);
D5 = nchoosek(v,6);
cov5=zeros(length(D5(:,1)),1);
cov=zeros(length(N2),1);
%This creates the matrix from the N2 N3 graph with only columns at the defenders and rows  N2. 
T=1:32;
for r=1:length(D5(:,1))
    DEFMAT=N2N3;
    DEFMAT(setdiff(T,N2),:)=[];
    DEFMAT(:,setdiff(T,D5(r,:)))=[]; 
    cov=zeros(length(N2),1);
        for s=1:length(N2)
        cov(s,1)=sum(DEFMAT(s,:));
        end
    cov5(r,1)=all(cov);
end
nnz(cov5)
%% Translation

