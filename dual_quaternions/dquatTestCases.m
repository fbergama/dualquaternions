classdef dquatTestCases < matlab.unittest.TestCase
    %DQUATTEST 
    
    properties
    end
    
    methods (Test)
        
        function testExpILogI( testCase )
            q = dquat_normalize(rand(8,1));
            testCase.verifyLessThan( sum( abs(dquat_ExpI(dquat_LogI(q)) - q)),1E-8 );
        end
        
        function testExpI2LogI2( testCase )
            q = dquat_normalize(rand(8,1));
            testCase.verifyLessThan( sum( abs(dquat_ExpI2(dquat_LogI2(q)) - q)),1E-8 );
        end
        
        
        function testExppLogp( testCase )
            q = dquat_normalize(rand(8,1));
            p = dquat_normalize(rand(8,1));
            testCase.verifyLessThan( sum( abs(dquat_Exp(p,dquat_Log(p,q)) - q)),1E-8 );
        end
            
    end
    
end

