classdef fluidDrivenAseismicSlip

    methods(Static)

        %% 2D problem
        
        function g = evalG2(xi,lambda,epsilon)
            a = sqrt(pi)*lambda.*exp(lambda.^2);
            erf1 = erf(lambda); erf2 = erf(lambda.*xi);
            g = (1+a.*(erf1-(1+epsilon).*erf2))./(1+a.*erf1);
            g(xi>1) = 0;
        end
        
        function T = evalT2(lambda,epsilon)
             %T = 2/pi*integral(@(xi) fluidDrivenAseismicSlip.evalG2(xi,lambda,epsilon)./sqrt(1-xi.^2),0,1);
             T = 2/pi*integral(@(z) fluidDrivenAseismicSlip.evalG2(cos(z),lambda,epsilon),0,pi/2);
        end
        
        function g = evalG2_BV(xi,lambda)
            g = erfc(lambda.*xi);
        end
        
        function T = evalT2_BV(lambda,symm)
            %T = 2/pi*integral(@(xi) fluidDrivenAseismicSlip.evalG2_BV(xi,lambda)./sqrt(1-xi.^2),0,1);
            T = 2/pi*integral(@(z) fluidDrivenAseismicSlip.evalG2_BV(cos(z),lambda),0,pi/2);
        end
        
        %% 3D problem
        
        function g = evalG3(xi,lambda,epsilon)
            g = expint(lambda.^2.*xi.^2)-expint(lambda.^2)+exp(-lambda.^2)./lambda.^2-epsilon;
            g(xi>1) = 0;
        end
        
        function T = evalT3(lambda,epsilon)
            %T = integral(@(xi) fluidDrivenAseismicSlip.evalG3(xi,lambda,epsilon).*xi./sqrt(1-xi.^2),0,1);
            T = integral(@(z) fluidDrivenAseismicSlip.evalG3(cos(z),lambda,epsilon).*cos(z),0,pi/2);
        end
        
        function g = evalG3_SL(xi,lambda)
            g = expint(lambda.^2.*xi.^2);
        end
        
        function T = evalT3_SL(lambda)
            %T = integral(@(xi) fluidDrivenAseismicSlip.evalG3_SL(xi,lambda).*xi./sqrt(1-xi.^2),0,1);
            T = integral(@(z) fluidDrivenAseismicSlip.evalG3_SL(cos(z),lambda).*cos(z),0,pi/2);
        end
        
    end

end