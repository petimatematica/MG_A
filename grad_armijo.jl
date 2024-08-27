using LinearAlgebra
using Printf

function descentgrad(x,  # Initial guess
                     f :: Function ,  # Objective function 
                     g :: Function ,  # Gradient of objective function
                     η :: Float64,  # Linesearch parameter
                     epsilon :: Float64, # Tolerance for gradient
                     maxiter :: Int64,  # Maximum of iterations 
                     stpmin :: Float64  # Minimum accetable step lenght 
                     ) 
   

                     # Iteration number initialization
                     iter = 0
                     stp = 0.0
                     while true

                        # Calculation of f at x
                        f_x = f(x)

                        # Calculation of gradient of function at x
                        gradf_x = g(x)

                        # Norm Calculation
                        norm_gradf_x = norm(gradf_x,2)

                        # Print info
                        @printf("%5d %20.15e %20.15e %20.15e\n",iter,norm_gradf_x,f_x,stp)

                        # Check for convergence
                        if norm_gradf_x < epsilon
                            return x,0
                        end

                        # Increment of iterations
                        iter += 1

                        # Check for maximum of iterations
                        if iter > maxiter
                            return x,1
                        end

                        # Armijo Linesearch
                        stp,ls_error = armijo(x,f,f_x,gradf_x,η,stpmin)

                        if ls_error > 0
                            return x,2
                        end

                        # Update sequence
                        x = x - stp * gradf_x
                     end
end


function armijo(x,f :: Function, f_x :: Float64, g_f_x, η :: Float64, stpmin :: Float64)

    stp = 1.0

    GD = -norm(g_f_x,2)^2 * η

    while true
        q = x - stp * g_f_x
        f_q = f(q)
        if f_q > f_x + stp *  GD
            stp = 0.5 * stp
            if stp < stpmin
                return stp,1
            end
        else
            return stp,0
        end

    end



end