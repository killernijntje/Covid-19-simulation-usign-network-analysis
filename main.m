clc; clear; 
%rng(3); % set to 3 to reproduce results  

plot_fr =1;
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% initial parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%MC repetions
reps=10;

%number of nodes
n=1000; 

% time limit
T=120;

% alternative model options

social_distancing =1;% 1 = no social distancing
lift_measures = 1; % 1 = no lift measures
lift_threshold = 0.1; % lift is treshold is lower than this 
central_hub = 1;% 1 = no hub
events = 1; % 1 = no events
event_limit = 100; % max no people in 1 event
graph_gif = 1; % 1 = no gif being created
herd =0; 
% SEIR parameters
    
%avg incubation time
avg_incub = 5.2;

%avg symptomatic time
avg_sympton=14; 
    
%basic reproductive number
R_0 = 2.2; 
%infection probability
eps= 1/avg_incub; 
gamma= 1/avg_sympton; 


%social distancing parameter
c1 = 0.01;
c2 = 0.05;
c3 = 0.1;



%counts
count_S = zeros(T,1,reps);
count_E = zeros(T,1,reps);
count_I = zeros(T,1,reps);
count_R = zeros(T,1,reps);

t_S1=zeros(reps,1);
t_S2=zeros(reps,1);
t_S3=zeros(reps,1);
lift_t=zeros(reps,1);


%social distancing checks
social=social_distancing;
removed_edges = [] ;

for l=1:reps
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% creating the network %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    q=1; 
    G = Waxman(n, q, round(n/n*10));
    avg_d = mean(degree(G));
    beta = R_0*(eps+gamma)/avg_d;
    G_plot=G;
    
    if social==0
        social_distancing=0;
    end
    check_s3=0;
    
    
    s1=0;
    s2=0;
    s3=0;
    %%%%%%%%%%%%%%%%%%%%%%%
    %%% SEIR SIMULATION %%%
    %%%%%%%%%%%%%%%%%%%%%%%
    S=ones(n,1);
    E=zeros(n,2);
    I=zeros(n,2);
    R=zeros(n,1);
    
    if herd==0
    vacinated = randperm(0.8*n); 
    S(vacinated)=0; 
    R(vacinated)=1;
    end
    %initial infected 
    int=randi([1,n]);
    S(int,1)=0;
    E(int,1)=1;
    E(int,2)=poissrnd(avg_incub) ;
    
    if graph_gif==0
        S_plot=zeros(T/plot_fr+1,n);
        E_plot=zeros(T/plot_fr+1,n);
        I_plot=zeros(T/plot_fr+1,n);
        R_plot=zeros(T/plot_fr+1,n);
        S_plot(1,:)=[G.Nodes.index(S==1); zeros(n-length(G.Nodes.index(S==1)),1)];
        E_plot(1,:)=[G.Nodes.index(E(:,1)==1) ;zeros(n-length(G.Nodes.index(E(:,1)==1)),1)];
    end
    for t=1:T
        if graph_gif==0
            if mod(t,plot_fr) == 0
                S_plot(t/plot_fr+1,:)=[G.Nodes.index(S==1); zeros(n-length(G.Nodes.index(S==1)),1)];
                E_plot(t/plot_fr+1,:)=[G.Nodes.index(E(:,1)==1) ;zeros(n-length(G.Nodes.index(E(:,1)==1)),1)];
                I_plot(t/plot_fr+1,:)=[G.Nodes.index(I(:,1)==1) ;zeros(n-length(G.Nodes.index(I(:,1)==1)),1)];
                R_plot(t/plot_fr+1,:)=[G.Nodes.index(R==1) ; zeros(n-length(G.Nodes.index(R==1)),1)];
            end 
        end


        [S,E,I,R] = SEIR(G,beta,S,E,I,R,avg_sympton,avg_incub);

        %check if central hub is 
        
        if central_hub == 0 
           cen_hub_ind = G.Nodes.inGraph(G.Nodes.inGraph>0);
           in_hub = randsample(cen_hub_ind,round(0.05*length(cen_hub_ind)));
           C = Waxman(length(in_hub),0.5,0);
           C.Nodes.hubindex(:)=in_hub;

           C_S=logical(S(in_hub)==1);
           C_E=ones(length(in_hub),2);
           C_E(:,1)=logical(E(in_hub)==1);
           C_I=ones(length(in_hub),2);
           C_I(:,1)=logical(I(in_hub)==1);
           C_R=logical(R(in_hub)==1);
           P_hub=R_0*(eps+gamma);%/mean(degree(C));
           [C_S_2,C_E_2,C_I_2,C_R_2] = SEIR(C,P_hub,C_S,C_E,C_I,C_R,avg_sympton,avg_incub);

           %Changing status after visit 
           S(in_hub) = S(in_hub) - (C_E_2(:,1)-C_E(:,1));
           E(in_hub,:) = E(in_hub,:) + C_E_2.*(C_E_2(:,1)-C_E(:,1));
        end 

    % check if events are enabled
        if events==0
            no_events = randi([1 10]);
            poss_to_attend = G.Nodes.inGraph(G.Nodes.inGraph>0);
            for i=1:no_events
                no_attendees=min(ceil(lognrnd(1,2))+1,event_limit);
                attending = randsample(poss_to_attend,min(no_attendees,length(poss_to_attend)));
                H = Waxman(length(attending),0.5,0);
                H.Nodes.eventindex(:)=attending;
                
                H_S=logical(S(attending)==1);
                H_E=ones(length(attending),2);
                H_E(:,1)=logical(E(attending)==1);
                H_I=ones(length(attending),2);
                H_I(:,1)=logical(I(attending)==1);
                H_R=logical(R(attending)==1);
                P_event=R_0*(eps+gamma);%/mean(degree(H));

                [H_S_2,H_E_2,H_I_2,H_R_2] = SEIR(H,P_event,H_S,H_E,H_I,H_R,avg_sympton,avg_incub);

                %Changing status after visit 
                S(attending) = S(attending) - (H_E_2(:,1)-H_E(:,1));
                E(attending,:) = E(attending,:) + H_E_2.*(H_E_2(:,1)-H_E(:,1));

                %allowing each node to visit only 1 event per period
                poss_to_attend(attending)=0;
                poss_to_attend(poss_to_attend==0)=[];
            end
        end


        E(E(:,2)>0,2)= E(E(:,2)>0,2)-1;% time until symptatic period reduced by one 
        I(I(:,2)>0,2)= I(I(:,2)>0,2)-1; % time until recory reduced by one 

        count_S(t,1,l) = sum(S);
        count_E(t,1,l) = sum(E(:,1));
        count_I(t,1,l) = sum(I(:,1));
        count_R(t,1,l) = sum(R);

        
        %lift measures
        if count_I(t)/n<lift_threshold && s3>6 && lift_measures == 0
            event_limit = 100; 
            social_distancing = 1; 
            EdgeTable = table(removed_edges, ones(length(removed_edges),1),...
                       'VariableNames',{'EndNodes','Weight'});
            G = addedge(G,EdgeTable);
            lift_t(l)=t;
        end
        
        
        
        
        %check if social distancing is applied
        if count_I(t)/n>c1 && social_distancing==0
            I_ind=find(I(:,1)==1);
            Isolate_s1 = randsample(I_ind,ceil(0.2*length(I_ind)));
            edges_s1_index=[];
            for j=1:length(Isolate_s1)
                edges_s1_index=[edges_s1_index ; ...
                                find(G.Edges.EndNodes(:,1)==Isolate_s1(j));...
                                find(G.Edges.EndNodes(:,2)==Isolate_s1(j))];
            end 
            removed_edges = [ removed_edges ; G.Edges.EndNodes(edges_s1_index,:)];
            G=rmedge(G,edges_s1_index);
            G.Nodes.inGraph(Isolate_s1)=0;
            s1=s1+1;
            if s1==1
                t_S1(l)=t; 
            end
        end

        
        
        if count_I(t)/n>c2 && social_distancing==0
            I_ind=find(I(:,1)==1); %
            select_I=G.Nodes.inGraph(I_ind);
            select_I(select_I>0)=1;
            I_ind=I_ind.*select_I; 
            I_ind(I_ind==0)= [];
            Isolate_s2 = randsample(I_ind,ceil(min(0.3+0.1*s2,0.7)*length(I_ind)));
            edges_s2_index=[];
            for j=1:length(Isolate_s2)
                edges_s2_index=[edges_s2_index ; ...
                                find(G.Edges.EndNodes(:,1)==Isolate_s2(j));...
                                find(G.Edges.EndNodes(:,2)==Isolate_s2(j))];
            end 
            removed_edges = [ removed_edges ; G.Edges.EndNodes(edges_s2_index,:)];
            G=rmedge(G, edges_s2_index);
            G.Nodes.inGraph(Isolate_s2)=0;
            s2=s2+1;
            if s2==1
              t_S2(l)=t; 
            end
            event_limit=20;
        end
% 
%         if count_I(t)/n>c3 && social_distancing==0
%             I_ind=find(S+E(:,1)+I(:,1)==1);
%             Isolate_s3 = randsample(I_ind,ceil((min(0.6+0.05*s3, 0.9))*length(I_ind)));
%             edges_s3_index=[];
%             for j=1:length(Isolate_s3)
%                 edges_s3_index=[edges_s3_index ; ...
%                                 find(G.Edges.EndNodes(:,1)==Isolate_s3(j));...
%                                 find(G.Edges.EndNodes(:,2)==Isolate_s3(j))];
%             end 
%             removed_edges = [ removed_edges ; G.Edges.EndNodes(edges_s3_index,:)];
%             G=rmedge(G,edges_s3_index);
%             G.Nodes.inGraph(Isolate_s3)=0;
%             s3=s3+1;
%             if s3==1
%               t_S3(l)=t; 
%             end
%             event_limit=2;
%             check_s3=1;
%         end
%         
        

        
        
    end

end 
%%


%%%%%%%%%%%%%%%%%%%%%
%%% plotting SEIR %%%
%%%%%%%%%%%%%%%%%%%%%


figure(2);
for l=1:reps
plot([0.4;count_S(:,1,l)./n], "LineWidth",1,"color",'#e8FFe8')
hold on
plot([0;count_E(:,1,l)./n],"LineWidth",1,"color",'#e8e8FF')
plot([0 ;count_I(:,1,l)./n],"LineWidth",1,"color",'#FFe8e8')
plot([0.6;count_R(:,1,l)./n],"LineWidth",1, "color",'#e8e8e8')
end 


plot([0.4;mean(count_S,[3])./n], "LineWidth",1,"color",'#00FF00')
% if social==0
%     if s1 > 0
%     xline(mean(t_S1), "m:");
%     end
%     if s2 > 0
%     xline(mean(t_S2), "m:");
%     end
%     if s3>0
%     xline(mean(t_S3), "m:");
%     end
%     xline(mean(lift_t), "r:");
% end
plot([0;mean(count_E,[3])./n],"LineWidth",1,"color",'#0000FF')
plot([0 ;mean(count_I,[3])./n],"LineWidth",1,"color",'#FF0000')
plot([0.6;mean(count_R,[3])./n],"LineWidth",1, "color",'#000000')


legend("S(t)","E(t)","I(t)","R(t)")
xlabel("Time (days)")
ylabel("fraction of population")

axis([0 T 0 1])
hold off


%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plotting graph transition %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if graph_gif==0
    h = figure;
    axis tight manual % this ensures that getframe() returns a consistent size
    filename = 'testAnimated2.gif';
    for k=1:T/plot_fr+1
        H = plot(G_plot); 
        highlight(H,S_plot(k,S_plot(k,:)>0),'NodeColor','g')
        highlight(H,E_plot(k,E_plot(k,:)>0),'NodeColor','b')
        highlight(H,I_plot(k,I_plot(k,:)>0),'NodeColor','r')
        highlight(H,R_plot(k,R_plot(k,:)>0),'NodeColor','k')
        drawnow

        frame = getframe(h); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256);   

          % Write to the GIF File 
          if k == 1 
              imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
          else 
              imwrite(imind,cm,filename,'gif','WriteMode','append'); 
          end 
    end
end 


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plotting graph transition %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if graph_gif==0
    h = figure;
    axis tight manual % this ensures that getframe() returns a consistent size
    filename = 'testAnimated_spatial.gif';
    X=G_plot.Nodes.xy(:,1);
    Y=G_plot.Nodes.xy(:,2);
    for k=1:T/plot_fr+1
        scatter(X,Y, 'filled')
        hold on
        scatter(X(S_plot(k,:)>0),Y(S_plot(k,:)>0),'filled','g')
        scatter(X(E_plot(k,:)>0),Y(E_plot(k,:)>0),'filled','b')
        scatter(X(I_plot(k,:)>0),Y(I_plot(k,:)>0),'filled','r')
        scatter(X(R_plot(k,:)>0),Y(R_plot(k,:)>0),'filled','k')
        hold off
        drawnow

        frame = getframe(h); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256);   

          % Write to the GIF File 
          if k == 1 
              imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
          else 
              imwrite(imind,cm,filename,'gif','WriteMode','append'); 
          end 
    end
end 