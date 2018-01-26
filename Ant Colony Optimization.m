function [bestAnt,pheromone,eta,objPlot,bottom]=ACO(arcMat,fleet,maxTime)
%ACO Generates a set of routes using ant colony optimization.
%   A set of routes is built by each ant in a colony, using a randomized
%   template as a cutoff limit and the beginning of a new route. Inputs are
%   an arcMat struct containing a 3D spare time matrix, a fleet struct with
%   robot parameters, and a maximum desired run time. Outputs the best
%   solution (routes) found, a matrix of peheromone values, a matrix of
%   heuristic values (eta), a record of the objective function over time,
%   and the lowest value in the initial spare time matrix (bottom).

ACOTimer=tic;
objPlot=zeros(10000,3); %store objective function values for plotting

%make all values in arc matrix >=0
bottom=min(arcMat.p(:));
if bottom<0
    arcMat.p=arcMat.p-bottom;
end

nNodes=size(arcMat.p,1);
nSensors=size(arcMat.p,2);
nRobots=size(arcMat.p,3);

%CHOOSE the following values as desired%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nAnts=100; %number of ants to use
Q=1/10000; %normalization constant
% R=1/10000;
pheromoneInitial=10e-6; %initial nonzero pheromone value
pie=1; %importance weight of pheromone value
sigma=1; %importance weight of heuristic value
rho=0.05; %pheromone evaporation rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%heuristic value of each element (arc) in the problem
%OPTION 1:
eta=Q*arcMat.p;
%OPTION 2:
% check=arcMat.p<-bottom; %all arcs with value <0 get low priority
% eta(check)=0.01;

%declare a pheromone matrix with uniform intial values
pheromone=pheromoneInitial*ones(nNodes,nSensors,nRobots);

%initialize an ant colony
emptyAnt=antRoutes;
emptyAnt.list=cell(1,fleet.size);
antColony=repmat(emptyAnt,nAnts,1);

%holder for the best set of routes found
bestAnt.minVal=-inf;
bestAnt.avgVal=-inf;

iter=0;
while toc(ACOTimer)<=maxTime
    
    iter=iter+1;
    bestMinAtIter=-inf;
    for a=1:nAnts %for each ant in the colony
        
        %set the robot start locations as the beginning of each route
        for k=1:nRobots
            antColony(a).list{k}=k+nSensors;
        end
        
        %generate a template of route loads for the current ant
        template=randomRouteTemplate(nSensors,fleet);
        
        %declare a binary list to represent visited sensor nodes
        visited=false(1,nSensors);
        
        for k=1:nRobots
            for l=1:template(k) %add number of sensors to route as per template
                
                i=antColony(a).list{k}(l); %set i to the end node in the current route
                
                %for each other node j in the problem, compute the
                %probability of selecting it as a successor to i, based on
                %a weighted combination of pheromone and heuristic value
                selectionProb=pheromone(i,:,k).^pie.*eta(i,:,k).^sigma;
                
                %set the probability of already visited nodes to zero
                selectionProb(visited)=0;
                
                selectionProb=selectionProb/sum(selectionProb);
                
                %use a roulette wheel to select the next node j, based on
                %the probabilities
                cumProb=cumsum(selectionProb);
                j=find(rand<=cumProb,1,'first');
                
                if j
                    antColony(a).list{k}(l+1)=j;
                    visited(j)=true;
                end
                
            end
        end
        
        %OPTIONAL: sort routes in each ant
%         for k=1:numRobots
%             antColony(a).list{k}(2:end)=sort(antColony(a).list{k}(2:end));
%         end

        antColony(a)=updateArcValuesAnt(antColony(a),arcMat);
        
        %update the best ant
        if antColony(a).minVal>bestAnt.minVal
            bestAnt=antColony(a);
            %OPTIONAL: Use the below if considering the max-average objective as well
%         elseif antColony(a).minVal==bestAnt.minVal
%             if antColony(a).avgVal>bestAnt.avgVal
%                 bestAnt=antColony(a);
%             end
        end
        
        if antColony(a).minVal>bestMinAtIter
            bestMinAtIter=antColony(a).minVal;
        end
        
    end
    
    %update pheromones
    %OPTION 1: Update pheromone from all ants
    for a=1:nAnts
        for k=1:nRobots
            arcs=[antColony(a).list{k}(1:end-1)' antColony(a).list{k}(2:end)'];
            pheromone_k=pheromone(:,:,k);
            phIndices=sub2ind(size(pheromone_k),arcs(:,1),arcs(:,2));
            pheromone_k(phIndices)=pheromone_k(phIndices)+Q*antColony(a).minVal/nAnts;%+R*antColony(a).avgVal/100;
            pheromone(:,:,k)=pheromone_k;
        end
    end

    %OPTION 2: Update pheromone from best ant only
%     for k=1:nRobots
%         arcs=[bestAnt.list{k}(1:end-1)' bestAnt.list{k}(2:end)'];
%         pheromone_k=pheromone(:,:,k);
%         phIndices=sub2ind(size(pheromone_k),arcs(:,1),arcs(:,2));
%         pheromone_k(phIndices)=pheromone_k(phIndices)+Q*bestAnt.minVal+R*bestAnt.avgVal;
%         pheromone(:,:,k)=pheromone_k;
%     end

    %evaporate pheromones
    pheromone=(1-rho)*pheromone;

    fprintf('ACO iteration %i best min: %f, best avg: %f, time: %f\n',iter,bestAnt.minVal,bestAnt.avgVal,toc(ACOTimer));
    objPlot(iter,1)=bestMinAtIter; %best min spare time from current iteration
    objPlot(iter,2)=bestAnt.minVal; %min spare time of best solution
%     objPlot(iter,3)=bestAnt.avgVal; %average spare time of best solution
    objPlot(iter,3)=toc(ACOTimer); %current run time

end

objPlot=objPlot(any(objPlot,2),:); %eliminate empty rows

end

