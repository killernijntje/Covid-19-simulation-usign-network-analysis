function G_A = Waxman(n, q, s)

%adjacency matrix (random network)
A = zeros(n);
location=rand(n,2);
distance = zeros(n);
for i=1:n
    for j=i+1:n
        distance(i,j)=sqrt( (location(i,1)- location(j,1))^2+(location(i,2)- location(j,2))^2 ) ;
        p_ij = q *exp(-s*distance(i,j)); 
        r=unifrnd(0,1);
        if r<= p_ij
            A(i,j)=1;
        else
            A(i,j)=0;
        end 
    end
end 
A = A - tril(A,-1) + triu(A,1)';
A = A - diag(diag(A))+ diag(zeros(n,1));
G_A = graph(A);
G_A.Nodes.xy = location;
G_A.Nodes.inGraph(:)=[1:n];
G_A.Nodes.index(:)=[1:n];
end 