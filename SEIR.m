function [S,E,I,R] = SEIR(G,p,S,E,I,R,avg_sympton,avg_incub)

% E -> I
I_ind = find(E(:,1)==1 & E(:,2)==0);
if isempty(I_ind)==0
    for i=1:length(I_ind)
        I(I_ind(i),:)=[1,poissrnd(avg_sympton)];
    end
end
E(I_ind,:)=0;


% I -> R

R(I(:,1)==1 & I(:,2)==0,:)=1;
I(I(:,1)==1 & I(:,2)==0,:)=0;

% S -> E
EI=E(:,1)+I(:,1);
EI_ind=find(EI==1);

%all nodes linked to a exposed or infected
pos_EI=[];
for i=1:length(EI_ind)
pos_EI=[pos_EI; G.Edges.EndNodes(G.Edges.EndNodes(:,1)==EI_ind(i),1:2); 
    G.Edges.EndNodes(G.Edges.EndNodes(:,2)==EI_ind(i),1:2)];
end 
%removing any exposed, infected or recovered nodes
EIR=E(:,1)+I(:,1)+R;
EIR_ind=find(EIR==1);
for i=1:length(EIR_ind)
    pos_EI(pos_EI==EIR_ind(i))=[];
end
%number of linkes betweens S and E/I
%unique posible EI
unique_pos_EI = unique(pos_EI);

infec_prob = p;

for i=1:length(unique_pos_EI)
    temp_node=unique_pos_EI(i);
    EI_near_nodes = sum(pos_EI(:) == temp_node);
    r=binornd(EI_near_nodes,infec_prob);
    if r>0
        S(temp_node,1)=0;
        E(temp_node,1)=1;
        E(temp_node,2)=poissrnd(avg_incub);
    end
end


end 
