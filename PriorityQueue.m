classdef PriorityQueue < handle
    
    % Priorigy queue data structure for the USC algorithm
    
    properties
        valueList;
        indexList;
        topindex;
        topvalue;
    end
    
    methods
        
        function thePQueue = PriorityQueue(n)
            thePQueue.valueList = sparse(1,n);
            thePQueue.indexList =[];
            thePQueue.topindex=0;
            thePQueue.topvalue=+inf;
        end
        
        function push(thePQueue, index, value)
            
            [index1,IA1] = setdiff(index, thePQueue.indexList);
            thePQueue.valueList(index1)=value(IA1);
            thePQueue.indexList = union(thePQueue.indexList,index1);
            
            
            [index2,IA2] = intersect(index, thePQueue.indexList);
            thePQueue.valueList(index2)=min(value(IA2),thePQueue.valueList(index2));
                        
            if min(value)<thePQueue.topvalue
                [thePQueue.topvalue,I]=min(value);
                thePQueue.topindex = index(I);
            end
        end
        
        
        function [minindex,minvalue] = pop(thePQueue)

                minindex = thePQueue.topindex;
                minvalue = thePQueue.topvalue;
                thePQueue.valueList(thePQueue.topindex)=0;
                thePQueue.indexList = setdiff(thePQueue.indexList,thePQueue.topindex);
                [~,J,V]=find(thePQueue.valueList);
                if ~isempty(J)
                    [thePQueue.topvalue,jindex] = min(V);
                    thePQueue.topindex=J(jindex);
                else
                    thePQueue.topindex=0;
                    thePQueue.topvalue=+inf;
                end
                
        end
        
        
        function flagIsEmpty = isEmpty(thePQueue)
            [~,J,V]=find(thePQueue.valueList);
            flagIsEmpty = isempty(J);
        end
        
    end
end