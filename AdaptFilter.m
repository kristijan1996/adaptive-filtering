classdef AdaptFilter
%%    
% Author: Kristijan Mitrovic
% Date: October 2018
% License: Open Source in every way
% 
% Description:
%     This class encapsulates basic least mean square (LMS)
%     algorithms used in adaptive filtering.
%     Every function is commented in detail in appropriate
%     place, this is just a brief description.
%     
%     During instantiation of the AdaptFilter object, user
%     passes two or three arguments to the constructor:
%     order of the filter, step size 'mu' used for weight
%     ajdustment and (optionaly) initial weights.
%     Something like this: 
%     >> filter = AdaptFilter(order, step, initW);
%     After that, he would call any of the available methods using
%     >> filter.method(signalToBeFiltered, desiredSignal);
%     The method will return the filtered signal, the error signal 
%     and the filter coeffitients.
%     
%     Few different methods have been implemented, all based on
%     basic LMS method.
%
%     IMPORTANT: Signals are passed and returned as row-vectors.
%%   
    properties (Access = private)
        L;      % Order of the filter, its Length
        mu;     % Step size 
        w;      % Filter coeffitients
    end
%%   
    methods
        
        function obj = AdaptFilter(order, stepSize, initialWeights)
            
            % AdaptFilter class constructor
            %
            % Call it like:
            % >> filter = AdaptFilter(N, mu, W0);
            % to instantiate an adaptive filter object of order
            % N with step size mu and initial weights W0.
            % Initial weights default to all zeros in case they
            % are not specified and the constructor is called
            % with just two arguments.
            
            if (nargin < 2 || nargin > 3)
                disp('Wrong number of input arguments!');
                return;
            elseif (nargin == 2)
                obj.L = order;
                obj.mu = stepSize;
                obj.w = zeros(1, obj.L);
            elseif (nargin == 3)
                obj.L = order;
                obj.mu = stepSize;
                obj.w = initialWeights;
            end               
                        
        end
%%       
        function [y, mse, w] = LMS(obj, x, d)
            
            % Calculates adaptive filter output y which
            % is obtained by filtering the input signal x, 
            % using signal d as a reference (desired) signal,
            % trying to minimize the difference e between the 
            % output and desired signal using gradient
            % descent algorithm.
            %
            % Call it using: 
            % >> [y, e, w] = filter.LMS(input, desired);
            % Where 'filter' is the instantiated adaptive
            % filter object.
            % It returns filtered signal y, error e between
            % desired and actual filter output and filter
            % coeffitients w.
             
            N = numel(x);
            y = zeros(1,N);
            e = zeros(1,N);
            
            % We need to pad input signal with L-1 zeros
            % because of Matlab indexing. From here on
            % we consider x_pad as actual input.
            % This needs to be done because the output
            % of the filter is calculated based on current
            % and past inputs and filter weights, and there 
            % are no past inputs before index 1 in the original
            % input sequence. Therefore, we need to "delay"
            % the input signal for L-1 samples.
            % Because of this, in places where we want to
            % get x(n), we will have to say x(n + L-1).
            
            x_pad = [zeros(1,obj.L-1), x];
            
            % Now for each input sample we calculate output,
            % error signal and update the filter weights
            for n = 1:N
                
                y(n) = obj.w * (x_pad(n+obj.L-1:-1:n))';
                
                e(n) = d(n) - y(n);
                mse(n) = 1/n*sum((e(1:n)).^2);
                
                obj.w = obj.w + obj.mu*e(n) * x_pad(n+obj.L-1:-1:n);
                
            end       
            mse = 1/N*sum(e);
            w = obj.w;
        end
%%      
        function [y, e, w] = signDataLMS(obj, x, d)
            
            % This is a modified version of the vanilla
            % LMS algorithm, which differentiates from
            % the original by sign(data) expression used when
            % updating filter weights.
            % For more info see the original LMS method.
            
            N = numel(x);
            y = zeros(1,N);
            e = zeros(1,N);
            
            x_pad = [zeros(1,obj.L-1), x];
            
            for n=1:N
               
                y(n) = obj.w * (x_pad(n+obj.L-1:-1:n))';
                
                e(n) = d(n) - y(n);
                
                obj.w = obj.w + obj.mu*e(n) * sign(x_pad(n+obj.L-1:-1:n));
                
            end       
            
            w = obj.w;
            
        end
%%        
        function [y, e, w] = signErrorLMS(obj, x, d)
            
            % This is a modified version of the vanilla
            % LMS algorithm, which differentiates from
            % the original by sign(error) expression used when
            % updating filter weights.
            % For more info see the original LMS method.
            
            N = numel(x);
            y = zeros(1,N);
            e = zeros(1,N);
            
            x_pad = [zeros(1,obj.L-1), x];
            
            for n=1:N
               
                y(n) = obj.w * (x_pad(n+obj.L-1:-1:n))';
                
                e(n) = d(n) - y(n);
                
                obj.w = obj.w + obj.mu* sign(e(n)) * x_pad(n+obj.L-1:-1:n);
                
            end       
            
            w = obj.w;
            
        end
%%        
        function [y, e, w] = signSignLMS(obj, x, d)
            
            % This is a modified version of the vanilla
            % LMS algorithm, which differentiates from
            % the original by sign(data)*sign(error) expression 
            % used when updating filter weights.
            % For more info see the original LMS method.

            N = numel(x);
            y = zeros(1,N);
            e = zeros(1,N);
            
            x_pad = [zeros(1,obj.L-1), x];
            
            for n=1:N
               
                y(n) = obj.w * (x_pad(n+obj.L-1:-1:n))';
                
                e(n) = d(n) - y(n);
                
                obj.w = obj.w + obj.mu* sign(e(n)) * sign(x_pad(n+obj.L-1:-1:n));
                
            end       
            
            w = obj.w;
            
        end
%%        
        function [y, e, w, mu] = NLMS(obj, x, d, alpha)
            % Normalized Least Mean Squares method
            %
            % This is a modified version of the LMS
            % algorithm which bypasses the issue of
            % having a constant step size in each step
            % of the algorithm.
            % The mu is calculated based on current 
            % energy of the input, which is proportional
            % to dot product of input by itself.
            % The parameter 0<alpha<2 is additionaly passed
            % as it helps tweak the convergence rate
            % of the algorithm. It defaults to 2 if
            % not specified, which results in max mu
            % for which the stability is guaranteed.
            
            if (nargin == 3)
                alpha = 2;
            end          

            N = numel(x);
            y = zeros(1,N);
            e = zeros(1,N);
            mu = 0.1*ones(1,N);
            
            x_pad = [zeros(1,obj.L-1), x];
            
            for n=1:N
               
                p = x_pad(n+obj.L-1:-1:n);
                
                y(n) = obj.w * p';
                
                e(n) = d(n) - y(n);
                
                if (alpha < 0 || alpha > 2)
                    mu(n) = 2/(obj.L*(p*p'));
                else
                    mu(n) = alpha/(obj.L*(p*p'));
                end
                
                obj.w = obj.w + mu(n) * e(n) * p;
                
            end       
            
            w = obj.w;
        end
        
    end
    
end