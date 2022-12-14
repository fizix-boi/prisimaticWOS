classdef Geometry
    properties
        COM
        DIM
        V
    end
    methods
        function obj = Geometry(com,dim,V)
            if nargin == 3
                obj.COM = com;
                obj.DIM = dim;
                obj.V   = V;
            end
        end
         function d = distance(obj,p)
             s = [0 0 0];
             for i = [1 2 3]
                 if p(i) < obj.COM(i)-obj.DIM(i)/2
                     s(i) = -1;
                 elseif  p(i) > obj.COM(i)+obj.DIM(i)/2
                     s(i) = 1;
                 else
                     s(i) = 0;
                 end
             end             
             pGeo = abs(s).*obj.COM + ( s.*obj.DIM./2 );
             pPoint = abs(s).*p;
             d = norm(pPoint-pGeo);
             if (length(find(s==0)) == 3)
                d = min(abs([p-(obj.COM+obj.DIM./2) p-(obj.COM-obj.DIM./2)]));
             end
         end
         function [] = PrintCorners(obj)
             fprintf('\n\n')
             fprintf('          [%1.2e,%1.2e,%1.2e]----------------[%1.2e,%1.2e,%1.2e]\n',[obj.COM(1)-obj.DIM(1)/2,obj.COM(2)+obj.DIM(2)/2,obj.COM(3)+obj.DIM(3)/2,obj.COM(1)+obj.DIM(1)/2,obj.COM(2)+obj.DIM(2)/2,obj.COM(3)+obj.DIM(3)/2]);
             fprintf('                    /|                                          /|\n');
             fprintf('                   / |                                         / |\n');
             fprintf('                  /  |                                        /  |\n');
             fprintf('                 /   |                                       /   |\n');
             fprintf('[%1.2e,%1.2e,%1.2e]----------------[%1.2e,%1.2e,%1.2e]\n',[obj.COM(1)-obj.DIM(1)/2,obj.COM(2)-obj.DIM(2)/2,obj.COM(3)+obj.DIM(3)/2,obj.COM(1)+obj.DIM(1)/2,obj.COM(2)-obj.DIM(2)/2,obj.COM(3)+obj.DIM(3)/2]);
             fprintf('                 |   |                                       |   |\n');
             fprintf('                 |   |                                       |   |\n');
             fprintf('                 |   |                                       |   |\n');
             fprintf('          [%1.2e,%1.2e,%1.2e]----------------[%1.2e,%1.2e,%1.2e]\n',[obj.COM(1)-obj.DIM(1)/2,obj.COM(2)+obj.DIM(2)/2,obj.COM(3)-obj.DIM(3)/2,obj.COM(1)+obj.DIM(1)/2,obj.COM(2)+obj.DIM(2)/2,obj.COM(3)-obj.DIM(3)/2]);
             fprintf('                 |  /                                        |  /\n');
             fprintf('                 | /                                         | / \n');
             fprintf('                 |/                                          |/  \n');
             fprintf('                 V                                           V   \n');
             fprintf('[%1.2e,%1.2e,%1.2e]----------------[%1.2e,%1.2e,%1.2e]\n',[obj.COM(1)-obj.DIM(1)/2,obj.COM(2)-obj.DIM(2)/2,obj.COM(3)-obj.DIM(3)/2,obj.COM(1)+obj.DIM(1)/2,obj.COM(2)-obj.DIM(2)/2,obj.COM(3)-obj.DIM(3)/2]);
             fprintf('\n')
         end
    end
end

