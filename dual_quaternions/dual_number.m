classdef dual_number
    %DUAL_NUMBER 
    
    properties
        a0=0
        ae=0
    end
    
    methods
        
        function dn = dual_number(nondual, dual)
            dn.a0 = nondual;
            dn.ae = dual;
        end
        
        function str = char(dn)
            
            sign='+';
            if dn.ae < 0
                sign='-';
            end            
            str = sprintf('%f %c %c%f',dn.a0,sign,char(400),abs(dn.ae));
        end
        
        
        function disp(obj)
           str = char(obj);
           disp(str);
        end 
        
        
        function dni = inv(dn)
            if dn.a0 == 0
                error('Unable to invert because non-dual part is zero');
            end
            dni = dual_number(1/dn.a0,-dn.ae/(dn.a0^2));
        end
        
        
        function dm = mtimes( dn1, dn2 )
            a0 = dn1.a0;
            ae = dn1.ae;
            b0 = dn2.a0;
            be = dn2.ae;
            
            dm = dual_number( a0*b0, a0*be + ae*b0 );
        end
        
        
        function dm = sqrt(dn)
            dm = dual_number( sqrt(dn.a0), dn.ae/(2*sqrt(dn.a0)));
        end
        
        
        function dr = sin(dn)
            dr = dual_number( sin(dn.a0), dn.ae*cos(dn.a0) );
        end
        
        
        function dr = cos(dn)
            dr = dual_number( cos(dn.a0), -dn.ae*sin(dn.a0) );
        end
        
        
        function dr = mrdivide( dn, num )
            dr = dual_number(dn.a0/num, dn.ae/num);
        end
        
    end
    
end

