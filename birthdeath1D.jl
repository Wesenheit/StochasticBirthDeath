using Plots,Random,Distributions

mutable struct Cell
    n::Int64 #number of elements
    l::Float64 #intensity of particle creation
    mu::Float64 #intensity of particle destruction(per element)
    t::Vector{Float64} #points of time
    N::Vector{Int64} #number of particles at given time
    Cell(a,b,c)=new(a,b,c,[],[])
end


function Evolve(system::Cell,num::Int64)
    push!(system.t,0)
    push!(system.N,system.n)
    for i in 1:num
        push!(system.t,rand(Exponential(system.l+system.n*system.mu),1)[1]+last(system.t))
        if rand()<= system.l/(system.l+system.mu*system.n)
            system.n+=1
        else
            system.n-=1
        end
        push!(system.N,system.n)
    end
end

function show(system::Cell)
    a=plot(system.t,system.N,
        title="Process of birth and death with \n λ=$(system.l) and μ=$(system.mu)",label="Number of particles",size=(1000,800))
    png(a,"out")
end


przyk=Cell(100,50,1)
Evolve(przyk,5000)
show(przyk)
