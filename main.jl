include("grad_armijo.jl")


# Função objetivo e seu gradiente
f(x) = 2*x[1]^2 - 1.05*x[1]^4 + (x[1]^6)/6 + x[1]*x[2] + x[2]^2
gradf(x) = [4*x[1] - 4.2*x[1]^3 + x[1]^5 + x[2], x[1] + 2*x[2]]

# Ponto inicial
x0 = [-2.0,3.0]

# Parâmetros para a função armijo
η = 1.e-4
stpmin = 1e-5
epsilon = 1.e-8
maxiter = 10000

x,error = descentgrad(x0,f,gradf,η,epsilon,maxiter,stpmin) 