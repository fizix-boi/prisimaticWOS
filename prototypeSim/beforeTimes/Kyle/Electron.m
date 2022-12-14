classdef Electron < handle
    
    properties
        home
        position
        charge
    end
    
    methods
        function obj = Electron(pos,char)
            obj.position = pos;
            obj.home = pos;
            obj.charge = char;
        end
        
        function step = walk(obj,R)
            theta = acos(2*rand()-1);
            phi = rand()*2*pi;

            step = [R*cos(phi)*sin(theta) R*sin(phi)*sin(theta) R*cos(theta)];
            obj.position = obj.position+step;
        end
        
        function [] = reset(obj)
            obj.position = obj.home;
        end
        
        function d = distance(obj)
            d = norm(obj.position - obj.home);
        end
        
        function [] = testWalk(obj,nRuns)
            figure(20)
            hold off
            for i = 1:nRuns
                step = obj.walk(1);
                scatter3(step(1),step(2),step(3));
                if mod(i,nRuns/100) == 0
                    fprintf('Run: %d \n',i)
                end
                hold on
            end
        end
    end
    
end

