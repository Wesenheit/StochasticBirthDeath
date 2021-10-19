using Plots,Random,Distributions

mutable struct Cell
    r::Int64 #number of m-RNA
    p::Int64 #number of proteins
    kr::Float64 #intensity of m-RNA creation
    gr::Float64 #intensity of m-RNA destruction(per element)
    kp::Float64 #intensity of protein creation
    gp::Float64 #intensity of protein decay(per element)
    t::Vector{Float64} #points of time
    R::Vector{Int64} #number of m-RNA at given time
    P::Vector{Int64} #number of proteins in time
    Cell(a,b,c,d,e,f)=new(a,b,c,d,e,f,[],[],[])
end

sample(items, weights) = items[findfirst(cumsum(weights) .> rand())]

function Evolve(system::Cell,num::Int64)
    push!(system.t,0)
    push!(system.R,system.r)
    push!(system.P,system.p)
    actions=["create_r","des_r","create_p","des_p"]
    for i in 1:num
        sum=system.kr+system.r*system.gr+system.p*system.gp+system.r*system.kp
        push!(system.t,rand(Exponential(sum),1)[1]+last(system.t))
        weights=[system.kr/sum,system.gr*system.r/sum,system.kp*system.r/sum,system.gp*system.p/sum]
        ac=sample(actions,weights)
        if ac=="create_r"
            system.r+=1
        elseif ac=="des_r"
            system.r-=1
        elseif ac=="create_p"
            system.p+=1
        else
            system.p-=1
        end
        push!(system.R,system.r)
        push!(system.P,system.p)
    end
end

function show(system::Cell)
    p1=plot(system.t,system.R,title="Number of m-RNA",label="Number of particles")
    p2=plot(system.t,system.P,title="Number of proteins",label="Number of partilces")
    pn=plot(p1,p2,layout=(2,1))
    png(pn,"plot")
end


przyk=Cell(100,0,100,2,10,5)
Evolve(przyk,10000)
show(przyk)
