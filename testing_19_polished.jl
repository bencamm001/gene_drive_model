#testing page #works on model 18
#################################################################
#set-up parameters
#################################################################



    #initialise some key variables
    #popsi is number of populations
    refresh(popi = 1000000, popsi = 2, indiv = 3, length_i = 50)

    #################################################################
    #population parameters
    #################################################################

    #number of generations to cycle over
    gens = 100

    #carrying capacity of area
    max_pop = 1000000

    #inbreeding coefficient
    Fis = fill(00, 1, 1, pops)

    #number of offspring per mating pair
    clutch = 10

    #################################################################
    #variables that change code base
    #################################################################
    migration = true
    dynamic = false
    mutation = false
    cleaning = false
    crash = false
    plotting = false
    histograms = false #plot conversion efficiency distributions
    sequence = "random" #"fasta" or "random"


    #################################################################
    #Selection variables
    #################################################################
    #whether or not selection pressure is present
    #in each population
    pressure = fill(1, 1, 1, pops) #1 is on 0 is off
    # pressure = reshape([1, 1], 1, 1, pops) #1 is on 0 is off

    #The percentage of a population exposed to selection pressure
    exposure = fill(1.0, 1, 1, pops)
    exposure[1] = 0.7

    if pops > 1
        exposure[2] = 0.2
    end


    #dominance of fitness cost
    dominance = fill(1.0, 3) #1 = dom, .5 = add, 0 = recess
    dominance[2] = 1
    dominance[3] = 1

    #random environmental affect on fitness
    env .= 0.0

    #fitness values based of sequence, arbitary full fitness sequence
    #relative fitness
    fitness = "set" #"seq" "relative"
    traits1 = [seq_traits(pop1, trait_distro); zeros(2, size(pop1, 2))]
    traits1[1,:] .= (traits1[1,:] ./ maximum(traits1[1,:]))

    #fitness cost of gene drive
    traits1[2,1] = -(0.45)
    traits1[2,2] = -0.0 #fitness cost of resistance to drive

    if fitness == "set"
        traits1[1,:] .= 1
    end

    if fitness == "relative" #untested#
        #trade of between resistance and fitness
        #as resistance increases, fitness decreases
        traits1[3,:] = (m_conver .- 1) ./ 5
        traits1[3,1] = 0
    end



    #################################################################
    #Sequence variables
    #################################################################
    #mutation rate of sequences
    mu_rate = 0

    #distribution of mutations
    gamma .= 1

    #################################################################
    #Conversion variables
    #################################################################

    conversion_distribution = "set" #"set" "bipolar" "beta" "normal" "sequence"

    if conversion_distribution == "bipolar"
        percent_peak1 = 0.5
        peak1 = 0.7
        se1 = 0.2
        peak2 = 0.2
        se2 = 0.2
        maxi = 1
        mini = 0.05
    end

    if conversion_distribution == "set"
        conv_set = 0.8
        res = 0.02
    end

    if conversion_distribution == "beta"
        shape1 = 3
        shape2 = 0.6
    end

    if conversion_distribution == "normal"
        nmean = 0.8
        nstd = 0.1
        maxi = 1
        mini = 0.01
    end

    if conversion_distribution == "sequence"
        grna = "ATGC"
        grna = split(grna, "")
        #affect of gRNA interaction
        priming = 1

    end


    if conversion_distribution == "bipolar"
        #bimodal conversion distribution
        #binomial to say where each is sampled from
        binom1 = rand(Binomial(1, percent_peak1), size(freq1,1))'
        binom2 = abs.(binom1 .- 1)

        norma = Normal(peak1, se1)
        m_conver1 = rand(norma, 1, size(pop1, 2)) .* binom1

        normb = Normal(peak2, se2)
        m_conver2 = rand(normb, 1, size(pop1, 2)) .* binom2

        m_conver = m_conver1 .+ m_conver2

        m_conver[m_conver .< 0] .= mini
        m_conver[m_conver .> 1] .= maxi
    end

    if conversion_distribution == "set"
        m_conver = fill(conv_set, 1, size(pop1, 2))
        m_conver[2] = res
    end

    if conversion_distribution == "beta"
        # artificial m_conver values
        conv_beta = Beta(shape1, shape2)
        m_conver = rand(conv_beta, 1, size(pop1, 2))
    end

    if conversion_distribution == "normal"
        # artificial m_conver values
        conv_norm = Normal(nmean, nstd)
        m_conver = rand(conv_norm, 1, size(pop1, 2))
        m_conver[m_conver .< 0] .= mini
        m_conver[m_conver .> 1] .= maxi
    end

    if conversion_distribution == "sequence"
        m_conver = m_conver_scan(pop1, drivers1, priming, grna)
    end

    if histograms == true
        histogram(m_conver', bins = 50, label = "", title = "initial distribution", xaxis = "conversion efficiency", yaxis = "frequency")
    end


    #################################################################
    #Population variation
    #################################################################

    #only works for 2 populations
    pop_variation = "none" #"subtle" "none"

    #the conversion efficiency at which the populations are separated
    cutoff = 0.5

    #the severity of the difference in populations
    magnitude = 10

    if pops == 2

        if pop_variation == "strict"
            freq1[:,:,1] .= freq1[:,:,1] .* transpose(m_conver .>  cutoff)
            freq1[:,:,2] .= freq1[:,:,2] .* transpose(m_conver .<= cutoff)
        end

        if pop_variation == "subtle"
            adjust1 = ones(size(freq1,1), 1)
            adjust1[m_conver' .>= cutoff] .*= magnitude
            adjust2 = ones(size(freq1,1), 1)
            adjust2[m_conver' .< cutoff] .*= magnitude

            freq1[:,:,1] .*= adjust1
            freq1[:,:,2] .*= adjust2
        end

        if pop_variation == "none"
            #do nothing
        end

    end

    if histograms == true
        ps = repeat([Plots.plot(1)], pops)
        ps[1] = scatter(m_conver',freq1[:,:,1], label = "", title = "Population 1", xaxis = "conversion efficiency", yaxis = "frequency")
        ps[2] = scatter(m_conver',freq1[:,:,2], label = "", title = "Population 2", xaxis = "conversion efficiency", yaxis = "frequency")
        plot(ps...)
    end


    #################################################################
    #Gene Drive starting frequency
    #################################################################

    #starting frequency of gene drive
    starting_freq = 0.001
    res_freq = 0.01
    val = sum(freq1[2:end,1,1]) + starting_freq
    freq1[1,1,1] = val * starting_freq
    freq1 ./= val

    # if pops > 1
    #     freq1[1,:,:] = [0.001, 0.0] #need array to be same size as pops
    # end

    freq1[1] = starting_freq
    freq1[2] = res_freq
    freq1[3] = 1 - freq1[1] - freq1[2]

    if pops > 1
        freq1[1,1,2] = 0
        freq1[2,1,2] = .3
        freq1[3,1,2] = 1 - freq1[1,1,2] - freq1[2,1,2]
    end



    #make all pop allele freqs sum to 1
    freq1 = mapslices(x -> x / sum(x), freq1, dims = [1,2])
    all1 = deepcopy(freq1)

    #################################################################
    #Migration rates
    #################################################################

    mig_style = "set" #"variable"

    if mig_style == "set"
        mig_rates .= 1e-4
    end

    if mig_style == "variable"
        mig_rates = reshape([0, .1, .1, 0], pops, pops) #vector needs to be pops squared length
    end


    #average trait value
    total_traits = Array{Float64}(undef, 0, pops)
    total_conv = Array{Float64}(undef, 0, pops)


    all1 = deepcopy(freq1)

# #change seed for simulation but keep starting scenario
# Random.seed!(50)

for q in 1:gens
    global freq1, pop1, all1, total_pop1, total_traits, sexes, sex_bias, pressure, dominance
    global popsize1, drivers1, traits1, m_conver, total_conv, total_sex, mig_rates, exposure


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


println(q)
end

# elapsed = round(time() - start, digits = 3)
# println(elapsed, " seconds")
###################################################################
#Plotting
###################################################################
#make empty array to store plots
pyplot()
ps = repeat([Plots.plot(1)], pops)

#make gene drives solid and others dashed
styles = fill(:dash, length(drivers1))
styles[vec(drivers1 .== 1)] .= :solid
styles = reshape(styles, 1, length(styles))

pop_names = ["Target", "Neighbour", "Souroukoudinga", "Ghana"]

for popu in 1:pops

    ps[popu] = Plots.plot(transpose(all1[:,:,popu]),
            #annotate = [(30, 0.3, text("Pressure = $(pressure[popu])", 8)),
            #            (30, 0.1, text("Env = $(env[popu])", 8))],
            legend = :topleft,
            xlab = "Generation \n",
            ylab = "Frequency",
            guidefontsize= 9,
            titlefontsize = 12,
            title = pop_names[popu],
            # left_margin = [-5mm 0mm],
            layout = 1,
            ylim = (0,1),
            linestyle = styles,
            width = 1.5,
            label = "")

        #     #plotting population size
        # Plots.plot!(total_pop1[:,:,popu] ./ max_pop,
        #     label = "",
        #     colour = "black",
        #     ls = :dash,
        #     width = 2)

        #     #plotting average fitness
        # Plots.plot!(total_traits[:, popu],
        # label = "",
        # colour = "black",
        # ls = :dot,
        # width = 2)

        #     #plotting average conversion
        # Plots.plot!(total_conv[:, popu],
        # label = "",
        # colour = "black",
        # ls = :dot,
        # width = 1.5)

        # Plots.vline!([200], color = "black", linewidth = 5, linestyle = :dot, label = "")
end




#plot all plots
Plots.plot(ps..., xlim = (0,gens), ylim = (0,1.05), reuse = false)







###################################################################
#Save file
###################################################################

# #save to file
# all1 = reshape(all1[1,:,:], 1, :, pops)
#
# @save "/home/student.unimelb.edu.au/bcamm/Documents/Uni/PhD/Modelling/Output/Test3/Test_run_c$(conversion)_m$(n)_t$(t)_rep$(m).jld2" all1 total_pop1
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
