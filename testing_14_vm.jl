#testing page
input = ARGS[6]
output = ARGS[7]
#output file path without terminal /

start = time()

include(input)


# Random.seed!(55)
# makes placeholder variables
refresh(popi = 1000000, popsi = 2, indiv = 2, length_i = 100)

index = parse(Int64, ARGS[1])
conv = parse(Float64, ARGS[2])
scost = parse(Float64, ARGS[3])
expo = parse(Float64, ARGS[4])
mig_input =parse(Float64, ARGS[5])




gens = 200
migration = true
max_pop = 1000000
dynamic = false
mutation = false
cleaning = false
crash = false

priming = 1


#distribution of mutations
gamma .= 1


clutch = 10
pressure[1] = 1
pressure[2] = 1
# pressure[3] = 0
# pressure[4] = 0
exposure = expo


#real Fis data

env .= 0.0
mu_rate = 0
#
# #put all the low conversion alleles in one population
# m_conver = m_conver_make(pop1, drivers1)

# # artificial m_conver values
# conv_norm = Normal(0.5, 0.2)
# m_conver = rand(conv_norm, 1, size(pop1, 2))
# m_conver[m_conver .< 0] .= 0.01
# m_conver[m_conver .> 1] .= 1

m_conver = ones(1, size(pop1, 2))
m_conver .= conv
# m_conver = m_conver_make(pop1, drivers1, priming, 1, 20, 31, 50)
# m_conver = m_conver_scan(pop1, drivers1, priming, grna)
# # histogram(transpose(m_conver))

# #strict haplotype position
# freq1[:,:,1] .= freq1[:,:,1] .* transpose(m_conver .>  .2)
# freq1[:,:,2] .= freq1[:,:,2] .* transpose(m_conver .<= .2)
# #

# #subtle haplotype position
# freq1[:,:,1] .= freq1[:,:,1] .+ (freq1[:,:,1] .* 10 .* transpose(m_conver .>  .5))
# freq1[:,:,2] .= freq1[:,:,2] .+ (freq1[:,:,2] .* 10 .* transpose(m_conver .<=  .5))


# freq1 = mapslices(x -> x / sum(x), freq1, dims = [1,2])

# freq1[freq1 .== 0] .= 0.00001
freq1 .= 1
freq1[1,1,1] = 0.001
freq1[1,1,2] = 0.00

# freq1[1,1,2] = 0.001
# freq1[4,1,2] = 0.00
# freq1[3,1,2] = 0.00

# freq1[1,1,3] = 0.001
# freq1[4,1,3] = 0.00
# freq1[2,1,3] = 0.00

# freq1[2,1,4] = 0.0
# freq1[3,1,4] = 0.0


#make all pop allele freqs sum to 1
freq1 = mapslices(x -> x / sum(x), freq1, dims = [1,2])
all1 = deepcopy(freq1)

# traits1[1,:] .= 0.45
# pop1[:,1] .= 1
traits1 = [seq_traits(pop1, trait_distro); zeros(2, size(pop1, 2))]
traits1[1,:] .= (traits1[1,:] ./ maximum(traits1[1,:])) ./ 2

traits1[1,:] .= 0.5
traits1[2,1] = -scost
# traits1[3,:] .= -0.4


# traits1[1,2] = 0.5
# m_conver[1,2] = 0.7
#
# traits1[1,3] = 0.4
# m_conver[1,3] = 0.01


mig_rates .= mig_input
# mig_rates[:,4] .= 0

# #JSON input files
# inputs = JSON.parse(String(read("input_markup.json")))
# inputs[

#average trait value
total_traits = Array{Float64}(undef, 0, pops)
total_conv = Array{Float64}(undef, 0, pops)


all1 = deepcopy(freq1)

# traits1[1,2] = 0.4
# m_conver[2] = .01
# freq1[2,:,:] .= 0
# freq1[2,:,4] .= 0.5
# freq1 = mapslices(x -> x / sum(x), freq1, dims = [1,2])

#
# #change seed for simulation but keep starting scenario
# Random.seed!(50)

for q in 1:gens
    global freq1, pop1, all1, total_pop1, total_traits, sexes, sex_bias
    global popsize1, drivers1, traits1, m_conver, total_conv, total_sex


    #1 generation of gene drive
    cycle1 = gene_drive(freq1, pop1, popsize1, drivers1, traits1, m_conver, Fis, sexes)
    #store new values for next cyle
    freq1 = cycle1[1]
    pop1 = cycle1[2]
    popsize1 = cycle1[3]
    drivers1 = cycle1[4]
    traits1 = cycle1[5]
    m_conver = cycle1[6]


if migration == true
    #migration
    #simulates how much of each allele travels to each pop
    migrants = repeat([freq1[:,:,1]],pops, pops)
    for k in 1:pops
        for j in 1:pops

            migrants[j,k] = freq1[:,:,j] .* mig_rates[j,k] .* popsize1[j]

        end
    end

    # #adjust each pops frequencies based on migration
    # #each row of migrants represents all the leavers
    # #each column represents all the arrivers
    for popu in 1:pops

        freq1[:,:,popu] = freq1[:,:,popu] .+ sum(hcat(migrants[:,popu]...) ./ repeat(fill(popsize1[popu][1], 1, size(popsize1, 3)), outer = size(pop1, 2)), dims = 2)

    end

    #normalise freq1
    freq1 = mapslices(x -> x / sum(x), freq1, dims = [1,2])
end



    #store new freq in all for plotting purposes
    all1 = combine_all1(all1, freq1, pops)


    #store total_pop
    total_pop1 = vcat(total_pop1, popsize1)
    total_sex = vcat(total_sex, sexes)


    #population cleaning
    #if allele is present for less than 10 generations,
    #remove it from:
    #pop1, freq1, drivers1, traits1, m_conver

        temp_conv = []
    for j in 1:pops
        push!(temp_conv, sum(m_conver .* repeat(transpose(freq1[:,:,pops]), size(m_conver, 1)) / size(m_conver, 1)))
    end
    total_conv = vcat(total_conv, transpose(temp_conv))


#cleaning dead alleles out
if cleaning == true
         #to remove
         to_keep = collect(1:size(all1,1))
         for k in 1:size(all1,1)

            if sum(all1[k, end, :] .< 1e-6) == pops &&
                sum(maximum(all1[k, :, :], dims = 1) .< 0.001) == pops

                deleteat!(to_keep, findall(to_keep .== k))

            end
        end



        #only keep the alleles specified in to_keep
        pop1 = pop1[:, to_keep]
        freq1 = freq1[to_keep, :, :]
        drivers1 = drivers1[: ,to_keep]
        traits1 = traits1[:, to_keep]
        m_conver = m_conver_make(pop1, drivers1, priming, 1, 20, 31, 50)
        all1 = all1[to_keep, :, :]
end



        temp_traits = []
    for j in 1:pops
        push!(temp_traits, 2 * sum(transpose(freq1[:,:,j]) .* sum(traits1, dims = 1)))
    end
    temp_traits[temp_traits .> 1] .= 1
    total_traits = vcat(total_traits, transpose(temp_traits))



    #dynamic bar plotting
    #plot data - make drives negative values
    #pdrives = transpose(map(x -> (-2 * x) + 1, drivers1))
    #pfreq = sort(freq1 .* pdrives, dims = 2)
if dynamic == true
    dat = transpose(freq1[:,1,:]) .* repeat(vec(popsize1), 1, size(freq1,1)) ./ max_pop
    nam = repeat(["Bana", "Pala"], outer = size(freq1,1))
    key = sort(unique(nam))
    val = collect(1:length(key)) .- 0.5
    dict = Dict(zip(key, val))
    ctg = repeat(1:size(freq1,1), inner = size(freq1, 3))
    p = groupedbar(nam, dat, bar_position = :stack,
        group = ctg,
        label = "", title = "Gene Drive Simulation \n Generation $(q)", linecolour = :match,
        ylim = (0,1.3), xlab = "Populations", ylab = "Frequency",
        ann = [ (dict["Bana"], 1.05, popsize1[1]),
                (dict["Pala"], 1.05, popsize1[2]),
                ])

    display(p)
    # savefig("frame$(q).png")
end


    # #dynamic pie charts
    # ps = repeat([Plots.plot(1)], pops)
    #
    # for popu in 1:pops
    #     ps[popu] = pie(freq1[:,:,popu])
    # end
    #
    # display(Plots.plot(ps...))


# println(q)
end

# elapsed = round(time() - start, digits = 3)
# println(elapsed, " seconds")
# ###################################################################
# #Plotting
# ###################################################################
# #make empty array to store plots
# ps = repeat([Plots.plot(1)], pops)
#
# #make gene drives solid and others dashed
# styles = fill(:dash, length(drivers1))
# styles[vec(drivers1 .== 1)] .= :solid
# styles = reshape(styles, 1, length(styles))
#
# pop_names = ["Bana", "Pala", "Souroukoudinga", "Ghana"]
#
# for popu in 1:pops
#
#     ps[popu] = Plots.plot(transpose(all1[:,:,popu]),
#             #annotate = [(30, 0.3, text("Pressure = $(pressure[popu])", 8)),
#             #            (30, 0.1, text("Env = $(env[popu])", 8))],
#             legend = :topleft,
#             xlab = "Generation \n",
#             ylab = "Frequency",
#             guidefontsize= 9,
#             titlefontsize = 12,
#             title = "$pop_names"[popu],
#             left_margin = [-5mm 0mm],
#             layout = 1,
#             ylim = (0,1),
#             linestyle = styles,
#             width = 1.5,
#             label = "")
#
#             #plotting population size
#         Plots.plot!(total_pop1[:,:,popu] ./ max_pop,
#             label = "",
#             colour = "black",
#             ls = :dash,
#             width = 2)
#
#             #plotting average fitness
#         Plots.plot!(total_traits[:, popu],
#         label = "",
#         colour = "black",
#         ls = :dot,
#         width = 2)
#
#         #     #plotting average conversion
#         # Plots.plot!(total_conv[:, popu],
#         # label = "",
#         # colour = "black",
#         # ls = :dot,
#         # width = 1.5)
#
#
#
# end
#
# #plot all plots
# Plots.plot(ps..., xlim = (0,gens), ylim = (0,1.05))
#
#





###################################################################
#Save file
###################################################################

# #save to file
# all1 = reshape(all1[1,:,:], 1, :, pops)
#
@save "$(output)/Test_run_c$(conv)_sc$(scost)_e$(expo)_m_$(mig_input)_rep$(index).jld2" all1 total_pop1
#
# end
# end
# end
# end
#
#



# #how to make unique identifiers for each sequence
# #works for very large bit strings
# mapslices(join, pop1, dims = 1)
# bytes2hex.(sha256.(mapslices(join, pop1, dims = 1)))
#







####
