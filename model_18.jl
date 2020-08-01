#making a model in julia 1.6.1
#
try
	#using these packages
	using Pkg, StatsBase, Compat, Random, Distributions, Plots, LinearAlgebra
	using Combinatorics, Random, LinearAlgebra, Base.Iterators, SHA
	using JLD2, CSV, RData, RCall, DataFrames, StatsPlots, Distributed, FastaIO
	using Measurements, LightXML, Statistics, PyPlot
catch
	using Pkg
	#add packages that are needed
	Pkg.add.(["RCall", "RData", "StatsBase", "Compat", "Random", "Distributions", "Plots", "LinearAlgebra"])
	Pkg.add.(["Random", "DataFrames", "CSV", "Combinatorics", "SHA", "JLD2", "StatsPlots", "Distributed", "FastaIO"])
	Pkg.add.(["Measurements", "LightXML", "Statistics", "PyPlot"])
	using Pkg, StatsBase, Compat, Random, Distributions, Plots, LinearAlgebra
	using Combinatorics, Random, LinearAlgebra, Base.Iterators, SHA
	using JLD2, CSV, RData, RCall, DataFrames, StatsPlots, Distributed, FastaIO
	using Measurements, LightXML, Statistics, PyPlot

end

#########################################################################################
#Function list
#########################################################################################

#############################
#initialise function
#############################
popi = 100
popsi = 4
migi = 0.0
conv = 0.8
indiv = 5
evni = 0
length_i = 10
refresh = function(;popi = 100, popsi = 4, migi = 0.0, conv = 0.8, indiv = 5, evni = 0, length_i = 10)

	#number of populations considered
	global pops
	pops = popsi

	#allele length
	global allele
	allele = length_i

	#starting alleles
	global init_allele
	init_allele = indiv

	#sequence of random haplotypes
	global pop1
	pop1 = convert(Array{Int8, 2}, rand(0:3, (allele, init_allele)));

	#frequency of haplotypes
	global freq1
	freq1 = round.(rand(init_allele, 1, pops), digits = 6);
		#normalise over 1
		freq1 = mapslices(x -> x / sum(x), freq1, dims = [1,2])

		#have to define this outside the function for some reason
	# #global frequency vector
	# global all1
	# all1 = deepcopy(freq1);

	#define genomic type
	global gene_type
	gene_type = repeat(["genomic"], size(pop1, 1));

	#define chromosome
	global chrm
	chrm = ones(size(pop1, 1));

	#define position on chromosome
	global pos
	pos = collect(1:size(pop1,2));

	#define population
	global pop
	pop = ones(size(pop1,2));

	#population size
	global popsize1
	popsize1 = reshape(repeat(popi:popi, pops), (1,1,pops));

	#maximum popsize
	global max_pop
	max_pop = 10000;

	#clutch size, affects maximum reproductive capacity of individual
	global clutch
	clutch = 100;

	#conversion effeciency of drive
	global conversion
	conversion = conv

	#define substitution matrix
	global sub_matrix
	#A = 0, T = 1, G = 2, C = 3, gap = 4
	sub_matrix = Dict(0 => [0.25, 0.25, 0.25, 0.25, 0], 1 => [0.25, 0.25, 0.25, 0.25, 0],
	 2 => [0.25, 0.25, 0.25, 0.25, 0], 3 => [0.25, 0.25, 0.25, 0.25, 0], 4 => [0, 0, 0, 0, 1]);

	#sex bias dictionary
	global sex_dict
	sex_dict = Dict(0 => 0.5, 1 => conversion, 2 => 1)


	#the mutation rate per generation
	global mu_rate
	mu_rate = 1e-5;

	#storage for all freqs in simulation
	global freq_all
	freq_all = freq1;

	#selection pressure 0 <= pressure <= 1
	global pressure
	pressure = zeros(1, 1, pops)

	#resistant haplotype that all others are compared to
	global res_haplo
	res_haplo = convert(Array{Int8, 2}, rand(0:3, (size(pop1, 1),1)))

	# #defines fitness of each haplotype against the reference
	# global haplo_fitness
	# haplo_fitness = mapslices(haplo_fit, pop1, dims = 1)

	#track total popsize
	global total_pop1
	total_pop1 = deepcopy(popsize1)

	#which haplotypes are drives
	global drivers1
	drivers1 = zeros(1, size(pop1, 2))
	drivers1[1] = 1

	#sex_bias
	global sexes
	sexes = fill(0.5, 1, 1, pops);

	#sex bias
	global sex_bias
	sex_bias = fill(0.5, 1, 1, pops);

	global total_sex
	total_sex = fill(0.5, 1, 1, pops);

	global trait_distro
	trait_distro = rand(Beta(2,100), size(pop1, 1))

	#multi fitness traits
	# global traits1
	# # traits1 = [rand(Normal(0.5, 0.1), 1, init_allele); zeros(2, init_allele)] ;
	# traits1 = [seq_traits(pop1, trait_distro); zeros(2, size(pop1, 2))]

	#migration matrix
	global mig_rates
	mig_rates = fill(migi, pops, pops)
	mig_rates[diagind(mig_rates)] .= 0

	#migration input
	global mig_input
	mig_input = migi

	#environmental factor
	global env
	env = zeros(1, 1, pops)

	# global mig_ref
	# mig_ref = (mig_a = mig_a, mig_b = mig_b, mig_c = mig_c, mig_d = mig_d, mig_e = mig_e)

	global exposure
	exposure = fill(1, 1, 1, pops)

	# global m_conver
	# m_conver = m_conver_make(pop1, drivers1)

	global ud_tox
	ud_tox = zeros(indiv, indiv)

	global Fis
	Fis = zeros(1,1,pops)

	global gamma
	gamma = ones(size(pop1, 1))

	global crash
	crash = false

	global dominance
	dominance = fill(1, 3)

end

##########################
#function end
##########################


#############################
#drift function - new frequecies with binomial
#############################

#define binomial distribution for change in frequency
#more complexity will be a unique distro for each haplotype give:
#popsize, clutch size

new_freq = function(x, popu)
	#need to for loop over popu
	#success of binomial should be eqaul to frequency
	#need to dot this function
	success = x

	distro = Binomial(2 * popsize1[:, :, popu][1], success)

	freq = rand(distro, 1)[1] / (2 * popsize1[:, :, popu][1])

end


##########################
#function end
##########################



#############################
#mutation function
#############################
mutate = function(haplo, position, sub_matrix, popu)

#haplo = the haplotypes that will mutate
#position = the position within the haplotype that will mutate
#sub_matrix = dictionary with mutation rates per base in order ATGC

#reset mu_haplo
	mu_haplo = Array{Int8}(undef, size(pop1, 1), 0)

	#input a haplotype array and a position array to mutate a base
	for i in 1:size(haplo, 2)

		#mutate single base in one haplo of input
		new_base = sample(0:3,
		Weights(sub_matrix[haplo[:,i][position[i]]]), 1)

		#check that the mutation is new, if so, add it
		if new_base != haplo[:,i][position[i]]
			temp_haplo = copy(haplo[:,i])
			temp_haplo[position[i]] = new_base[1]
			mu_haplo = hcat(mu_haplo, temp_haplo)
		end
	end

	return (mu_haplo)
end

##########################
#function end
##########################


##########################
#use mutate function
##########################

do_mutate = function(population, popsize, mu_rate)

	#number of mutations
	mu_number = trunc(Int, ceil(size(population, 1) * popsize[:, :, popu] * mu_rate))

	#column and row to be mutated in pop1
	mu_col = sample(1:size(population, 2), mu_number)
	mu_row = sample(1:size(population, 1), mu_number)
	#find haplotypes that mutated, copy them with mutated base
	#and add a new frequency of 1/2N

	#find the column where mutation occured, div divides without remainder
	#add copy of mutated haplotype to population array
	#get columns to mutate
	new_haplo = population[:, mu_col, popu]
	sites = mu_row

	#mutate the desired bases and return an array of new haplos
	mu_haplo = mutate(new_haplo, sites, sub_matrix)

end


#improved mutation module

new_mute = function(pop1, popsize1, mu_rate, sub_matrix, drivers1, gamma)

	#new alleles
	# global new_haplos
	new_haplos = Array{Int8}(undef, size(pop1,1), 0)
	new_freqs = Array{Float64}(undef, size(popsize1, 3), 0)
	new_drivers = Array{Float64}(undef, 1, 0)
	new_traits = Array{Float64}(undef, 3, 0)
	#mutations per population
	for popu in 1:pops

		#number of mutations in population
		mu_number = rand(Binomial(popsize1[popu], mu_rate), 1)[1]

		# #location of mutations - doesnt account for allele freqs
		# mu_locations = sample(1:prod(size(pop1)), mu_number)

		#alleles where mutations occur
		mu_alleles = sample(1:size(pop1, 2), Weights(vec(freq1[:,:,popu])), mu_number)

		#position of mutations along allele
		mu_position = sample(1:size(pop1, 1), Weights(gamma), mu_number)

		#combine position and allele
		location = hcat.(mu_position, mu_alleles)

		for i in location

			#copy the column allele and mutate
			new_allele = copy(pop1[:,i[2]])
			old_base = new_allele[i[1]]
			new_base = sample(0:4, Weights(sub_matrix[old_base]), 1)[1]

			if old_base != new_base
				new_allele[i[1]] = new_base

				#global new_haplos, new_freqs, new_drivers #needed if working inside function
				new_haplos = hcat(new_haplos, new_allele)

				#add new frequency but only to the right population
				temp_freq = zeros(size(popsize1,3))
				temp_freq[popu] = 1 / (2 * popsize1[popu])
				new_freqs = hcat(new_freqs, temp_freq)

				#add new haplos to drivers1 array
				new_drivers = hcat(new_drivers, drivers1[i[2]])

				# #add new trait values deviating from original
				# temp_traits = traits1[:, i[2]] .+ [rand(Normal(0, 0.05), 1); 0; 0]
				# new_traits = hcat(new_traits, temp_traits)

				#add new trait values based on sequence
				temp_traits = seq_traits(new_allele, trait_distro)[1]
				former_traits = copy(traits1[:, i[2]])
				former_traits[1] = temp_traits
				new_traits = hcat(new_traits, former_traits)

			end

		end

	end
	return(new_haplos, reshape(transpose(new_freqs), :, 1, size(popsize1,3)), new_drivers, new_traits)
end


#generate parameters for new sequences

##########################
#generate new trait values - same for all pops
##########################

new_trait = function(new_haplos, popsize1)

	#generate random trait values for new haplos
	new_traits = rand(Normal(.035, .01), 3, size(new_haplos, 2))
	#new_traits = fill(0.05, 3, size(new_haplos, 2), size(popsize1, 3))
end

seq_traits = function(pop1, trait_distro)

traits_1 = []

	for i in 1:size(pop1, 2)
		len = size(pop1, 1)
		push!(traits_1, sum((pop1[1:end,i] .== zeros(len, 1)) .* trait_distro))
	end
	return(transpose(traits_1))
end






##########################
#function end
##########################


##########################
#concatenate frequencies with preceding zeros
##########################

#old function
function cat_freq(freq_new, freq_old)

	#find number of new haplos compared to previous
	no_haplos = length(freq_new) - size(freq_old, 1)

	#fill an array of zeros for the number of generations and number of new haplos
	gap_fill = fill(0.0, no_haplos, size(freq_old, 2))

	#add to to the bottom of old frequency
	freq_old = vcat(freq_old, gap_fill)

	#cat the two freqs together
	freq_old = hcat(freq_old, freq_new)

end

#new function combine all1 and freq1

combine_all1 = function(all1, freq1, pops)

	#find how many new alleles there are
	newbies = size(freq1, 1) - size(all1, 1)

	#add zeros to all1 to accomodate new alleles
	all1_new = cat(dims = 1, all1, repeat(zeros(newbies, size(all1, 2)), outer = [1,1,pops]))

	all1_new = cat(dims = 2, all1_new, freq1)

	all1_new
end


##########################
#function end
##########################


##########################
#fitness module
##########################

haplo_fit = function(population)

	#input is population haplotypes, output is fitness of haplotypes
	#define fitness for each haplotype = likihood to survive
	#currently ratio of similarity to res_haplo is fitness
	#need to mapslices function

	fitness = sum(res_haplo .== population) / length(res_haplo)

end

##########################
#end
##########################

##########################
#selection diploids
##########################

#fitness of each pair of haploids
#selection is additive of each haplotype
selection_dips_add = function(fitness)

	#length of fitness
	len = length(fitness)

	#repeated fitness
	fit = repeat(fitness, len)

	#sum haplo fits for each combination of haplos
	dip_fit = fit + rotr90(fit)

	#fitness > 1 should be returned to 1
	dip_fit[dip_fit .> 1] .= 1

	return(LowerTriangular(dip_fit))

end


#fitness of each pair of haploids
#selection is dominant of each haplotype
selection_dips_dom = function(fitness)

	#length of fitness
	len = length(fitness)

	#repeated fitness
	fit = repeat(fitness, len)

	#cat the two arrays then take the max of each position
	temp = cat(fit, rotr90(fit), dims = 3)
	dip_fit = reshape(maximum(temp, dims = 3), (len, len))

	#fitness > 1 should be returned to 1
	dip_fit[dip_fit .> 1] .= 1

	return(LowerTriangular(dip_fit))

end

##########################
#multi-trait selection
##########################

# #how many alleles
# len = length(freq1)
#
# #index position of each diploid
# dips = collect.(collect(product(1:len, 1:len)))

#select each diploid and calculate fitness for many traits
#index the traits by the dips index and sum the columns for now
	trait_calc_add = function(dips, traits1, popu, dominance)
		#need to dot this function
		trait_score = Array{Float64}(undef, size(dips))
		for i in vec(dips)
			#base fitness
			base = mean(traits1[1, i])

			#fitness cost homos trait2
			if traits1[2, i][1] == traits1[2, i][2]
				cost1 = traits1[2, i][1]
			end

			#fitness cost heteros trait2
			if traits1[2, i][1] != traits1[2, i][2]
				cost1 = minimum(traits1[2, i]) * dominance[2]
			end

			#fitness cost homos trait3
			if traits1[3, i][1] == traits1[2, i][2]
				cost2 = traits1[3, i][1]
			end

			#fitness cost heteros trait3
			if traits1[3, i][1] != traits1[2, i][2]
				cost2 = minimum(traits1[3, i]) * dominance[3]
			end


			trait_score[i[1], i[2]] = sum([base, cost1, cost2])

			# trait_score[i[1], i[2]] = sum(traits1[:, i]) + (randn(1) * env[:,:,popu][1])[1]
		end
	trait_score
	end

	trait_calc_dom = function(dips, traits1, popu)
		#need to dot this function
		trait_score = Array{Float64}(undef, size(dips))
		for i in vec(dips)
		trait_score[i[1], i[2]] = maximum(sum(traits1[:, i], dims = 1)) + (randn(1) * env[:,:,popu][1])[1]
		end
	trait_score
	end



##########################
#end
##########################


#sample how many survive after selection
survival = function(dip_freq, dip_fit, popu, exposure)

	#survival(frequency of diploids, fitness of diploids)
	#sample how many expect to survive due to selection
	#binomial with trials = number of diploid, success = fitness
	#need to dot this function
	#survival is now how many die
	trials = trunc(Int, floor(exposure * dip_freq * popsize1[:,:,popu][1]))
	success = (1 - (dip_fit)) * pressure[:,:,popu][1]
	distro = Binomial(trials, success)
	surviving = rand(distro, 1)[1] / (popsize1[:,:,popu][1])

end


##########################
#make haploids from diploids
##########################

dip2hap = function(diploids)
	#how to return from diploid freqs to haploid freqs
	hap_freqs = []
	for i in 1:size(diploids, 1)

		#sum along the rows and columns and subtract the one counted twice
		new_freq = (sum(diploids[i,:]) / 2) + (sum(diploids[:,i]) / 2)
		 			+ diploids[i,i]

		if new_freq < 1e-20
			new_freq = 0
		end

		push!(hap_freqs, new_freq)

	end

	return(hap_freqs)

end


##########################
#function end
##########################


##########################
#find drivers
##########################

	driver_split = function(drivers)

		#separate each 1 in drivers into a unique array / new row
		len = length(drivers)
		base = zeros(Int8, len, len)

		#add the one for each driver in a new row
		split_drive = base .+ Diagonal(vec(drivers))

		#only keep the rows where drivers is 1
		split_drive = split_drive[Bool.(vec(drivers)), :]

	end



	driver_inter = function(drivers)
	#number of haplotypes

	drivers .+ transpose(drivers)

	end

##########################
#function end
##########################


##########################
#gene drive conversion - change frequencies
##########################

	converter = function(diploids, het_drivers)

		#frequency of hets should decrease and freq of
		#homo drive should increase

		heteros = map(x -> 1 - (x - 1)^2, het_drivers)
		homos = map(x -> x * (x - 1) * (x - 1.5), het_drivers)

		#find heterozygotes and record the change in their freqs
		change = diploids .* heteros .* conversion

		#add change to freq of homo drive
		diploids = diploids + (homos .* sum(change)) - change


		return(diploids)

	end



	converter1 = function(diploids, split_drive, m_conver, popu, drivers1)

		#get hetero/homo for gene drive
		het_drivers = LowerTriangular(split_drive .+ transpose(split_drive))
		heteros = map(x -> 1 - (x - 1)^2, het_drivers)

		position = findall(split_drive .==1)[1]
		drive_number = sum(drivers1[:, 1:position])
		conversions = transpose(m_conver[trunc(Int64, drive_number), :])

		#find opposite of drivers1 so can make
		#it that drivers cant convert drives
		not_drivers = (map(x -> -x + 1, drivers1))

		#make it relative conversion????
		conversions = conversions ./ maximum(conversions)
		conversions = conversions .* not_drivers #.* 2.5

		#need to remove the position of the homozygote gene drive
		homo_drive = findall(split_drive .== 0)
		conversions = conversions[:, homo_drive]

		targets = findall(heteros .== 1)
		convertees = diploids[targets]
		convert_rates = conversions
		conv_change = zeros(length(convertees))

		for i in 1:length(convertees)

			trials = trunc(Int64,  popsize1[popu])
			success = convertees[i] * convert_rates[i]
			distro = Binomial(trials, success)
			conv_change[i] = rand(distro, 1)[1] / popsize1[popu]
		end

		heteros[findall(heteros .== 1)] = heteros[findall(heteros .== 1)] .* conv_change
		return(heteros)
	end


converter2 = function(diploids, split_drive, m_conver, popu, drivers1)

	#get hetero/homo for gene drive
	het_drivers = LowerTriangular(split_drive .+ transpose(split_drive))
	heteros = map(x -> 1 - (x - 1)^2, het_drivers)
	targets = het_drivers .> 0

	#position of gene drive allele
	homo_pos = findall(split_drive .== 1)[1]
	drive_num = trunc(Int64, sum(drivers1[1:homo_pos]))

	#opposite of drivers1
	wilds = (-1 .* drivers1) .+ 1

	#conversion rate per allele
	#removing chance to convert other drives
	conversion = transpose(m_conver[drive_num, :]) .* wilds

	convertees = diploids .* targets
	convertees1 = convertees[targets]

	#conversion process
	conv_change = zeros(length(convertees1))
	for i in 1:length(convertees1)

		trials = trunc(Int64,  popsize1[popu] )
		success = conversion[i] * convertees1[i]
		distro = Binomial(trials, success)
		conv_change[i] = rand(distro, 1)[1] / popsize1[popu]
	end

	convertees[targets] .= conv_change

	return(convertees)

end






m_conver_make = function(pop1, drivers1, priming, g_start, g_end, t_start, t_end)

	#make empty array
	m_conver = Array{Float64}(undef, sum(drivers1 .== 1), size(pop1, 2))
	len = length(pop1[g_start:g_end,1])

	drivers = findall(drivers1 .== 1)
	for i in 1:sum(drivers1 .== 1)

		#only pick the drive alleles
		this_one = drivers[i][2]

		m_conver[i, :] = sum(pop1[g_start:g_end, this_one] .== pop1[t_start:t_end,:], dims = 1) ./ len .* priming

	end
	return(m_conver)
end


m_conver_scan = function(pop1, drivers1, priming, grna_seq)

	#make empty array
	m_conver = Array{Float64}(undef, sum(drivers1 .== 1), size(pop1, 2))
	len = length(grna_seq)

	drivers = findall(drivers1 .== 1)


	for i in 1:sum(drivers1 .== 1)

		#only pick the drive alleles
		this_one = drivers[i][2]

		grna = pop1[1:len, this_one]

		#iterate over sequence
		for k in 1:size(pop1, 2)

			mean_conv = []

				for j in 1:(size(pop1, 1) - len)

					push!(mean_conv, sum(grna .== pop1[j:(j + len - 1), k]) / len * priming)

				end

			m_conver[i, k] = maximum(mean_conv)

		end



	end
	return(m_conver)
end




	under_d = function(diploids, ud_tox)

		toxins = diploids .* ud_tox

		toxins = mapslices(x -> x ./ sum(x), toxins, dims = [1,2])

		toxins
	end






##########################
#function end
##########################


##########################
#sex bias conversion - male bias
#needs work!
##########################

#expected ratio of M/F offspring

	sex_conv = function(diploids, het_drivers)

		#diploids is the frequency of diploids
		#drivers is the homo/hetero state of drives

		#number of trials in binomial
		trials = trunc(Int, diploids * popsize)

		#success rate of binomial = ratio of M/F
		success = sex_dict[het_drivers]

		#Binomial distribution to see how many males born
		distro = Binomial(trials, success)

		sex = rand(distro, 1)[1] / trials



	end

##########################
#function end
##########################


##########################
#make diploids
##########################
make_dips = function(freq1, Fis)
	#need to for loop over pops
	#make a square matrix of freqs
	#make diploids with frequency of each diploid
	diploids = freq1 * transpose(freq1)

	#multiple everything except the diagonal by 2 to compensate for lower
	heteros = map(x -> -x + 1, Diagonal(ones(length(freq1))))
	homos = abs.(heteros .-1)

	#adjust for inbreeding
	diploids_Fis = Diagonal(vec(sum(diploids .* heteros .* (Fis), dims = 2)))
	diploids_homo = diploids .* homos .+ diploids_Fis
	diploids_hets = LowerTriangular(2 .* diploids .* heteros .* (1 - Fis))

	diploids = diploids_homo .+ diploids_hets
	diploids = diploids ./ sum(diploids)

end



##########################
#function end
##########################




##########################
#migration sampling
##########################

mig_sample = function(freq1, rate)


	#migration(frequency of haploids, popsize)
	#sample how many expect to survive due to selection
	#binomial with trials = number of diploid, success = fitness
	#need to dot this function
	trials = trunc(Int, round(from * fpop))
	success = rate
	distro = Binomial(trials, success)
	migrating = rand(distro, 1)[1] / fpop

end


##########################
#function end
##########################



##########################
#death rate per group
##########################
#returns the number of survivors

death = function(d_rate, pops)

    #function to see how many of each group die based on the
    #death rate ~ popsize ~ carrying capacity
    #need to dot this function

    #distribution for chance of survival (1 - death)
    ded = Binomial(pops, (1 - d_rate))

    survivors = rand(ded, 1)[1]

end


##########################
#gene drive total function
##########################

gene_drive = function(freq1, pop1, popsize1, drivers1, traits1, m_conver, Fis, sexes)

    # global freq1
    # global all
    # global pop1
    # global haplo_fitness
    # global total_pop
    # global popsize
    # global sexes
    # global drivers1
    # global traits
	# global m_conver

###################################################################
#mutation
###################################################################

if mutation == true

	#mutate sequences to make new alleles
	new_haplos, new_freqs, new_drivers, new_traits = new_mute(pop1, popsize1, mu_rate, sub_matrix, drivers1, gamma)

	#add sequences to pop1
	pop1 = hcat(pop1, new_haplos)

	#add new frequencies to freq1
	freq1 = cat(dims = 1, freq1, new_freqs)
	freq1 = mapslices(x -> x / sum(x), freq1, dims = [1,2])

	#adjust the drivers1 array
	drivers1 = hcat(drivers1, new_drivers)

    #add new traits to array - multi-trait fitness
    #traits1 = cat(dims = 2, traits1, new_trait(new_haplos, popsize1))
	# traits1 = hcat(traits1, new_traits)
	temp = seq_traits(pop1, trait_distro)
	traits1 = hcat(traits1, new_traits)
	traits1[1,:] = temp ./ maximum(temp) ./ 2


	# # #add new conversion efficiency (by similarity)
	# m_conver = m_conver_make(pop1, drivers1, priming, 1, 20, 31, 50)

	#scanning m_conver
	m_conver = m_conver_scan(pop1, drivers1, priming, grna)
end



###################################################################
#diploids/selection
###################################################################

    #make diploids with frequency of each diploid, sums to 1
    diploids = zeros(size(freq1,1), size(freq1,1), pops)
    for popu in 1:pops
        diploids[:,:,popu] = make_dips(freq1[:,:,popu], Fis[:,:,popu][1])
    end
#######################################################################
#######################################################################

###################################################################
#gene drive action - conversion
###################################################################
# drivers1 = zeros(size(freq1,1))
# drivers1[1] = 1
    #split the drivers array so each gene drive is solo
    #then iterate over each driver array
    #add the one for each driver in a new row
    split_drive = driver_split(drivers1)


#repeat for each row of split_drive the conversion process
for popu in 1:pops
	for i in randperm(size(split_drive,1))

        #convert het drives to homo drives
        diplo_change = LowerTriangular(converter2(diploids[:,:,popu], split_drive[i,:], m_conver, popu, drivers1))

		diploids[:,:,popu] = diploids[:,:,popu] .- diplo_change
        diploids[:,:,popu] = diploids[:,:,popu] .+ Diagonal(vec(sum(diplo_change) .* split_drive[i,:]))

		#make sure no values below 0 or above 1
		diploids[:,:,popu] = abs.((diploids[:,:,popu] .> 0) .* diploids[:,:,popu])
		diploids[:,:,popu] = diploids[:,:,popu] ./ sum(diploids[:,:,popu])

    end
end

###################################################################
#selection
###################################################################

    # #matrix of fitness for each diploid
    # fitness = selection_dips_add(haplo_fitness)

    #how many alleles
    len = length(freq1[:,:,1])

    #index position of each diploid
    dips = collect.(collect(product(1:len, 1:len)))

    #matrix of diploid fitness for multi-traits
	fitness = zeros(len, len, pops)
	for popu in 1:pops
		fitness[:,:,popu] = LowerTriangular(trait_calc_add(dips, traits1, popu, dominance))
	end

	fitness = fitness ./ maximum(fitness)


    # #set max fitness to 1 and min to 0
    # fitness[fitness .> 1] .= 1
	# fitness[fitness .< 0] .= 0

    #survival - see how many survive by fitness
	#survival is now how many die
    survivors = zeros(size(freq1,1), size(freq1,1), pops)
    for popu in 1:pops
        survivors[:,:,popu] = survival.(diploids[:,:,popu], fitness[:,:,popu], popu, exposure[:,:,popu])
    end

	survivors = diploids .- survivors
	survivors[survivors .< 0] .= 0


    #modify population size by how many died
	for popu in 1:pops
		popsize1[popu] = trunc.(Int64, popsize1[:,:,popu][1] * sum(survivors[:,:,popu]))
	end

	#issues arise if popsize = 0, make it equal to 1?
	popsize1[popsize1 .<= 0] .= 1

    #diploids frequency is survivors
    diploids = survivors

    #return diploids to sum = 1?
    diploids = mapslices(x -> x ./ sum(x), diploids, dims = [1,2])


###################################################################
#gene drive action - conversion sex bias
###################################################################

if crash == true
	for popu in 1:pops
	    #ratio of sexes
	    ratios = LowerTriangular(fill(0.5, size(diploids[:,:,popu])))

	    #change driver to be 1 (all males)
		het_drivers = drivers1 .+ transpose(drivers1)
		ratios[Diagonal(het_drivers .== 2)] .= 1

	    #convert sex ratio to average sex bias
	    sex_bias[:,:,popu] .= sum(diploids[:,:,popu] .* ratios)


	end

	#add sex bias to plotting vector
	sexes = vcat(sexes, sex_bias)
end


###################################################################
#Return to haploids for next generation
###################################################################

    #return the survivor freqs to haploid freqs
    for popu in 1:pops
        freq1[:,:,popu] = dip2hap(diploids[:,:,popu])
    end


    # #return new_freq to sum of 1
    # freq1 = freq1 ./ sum(freq1)

    #increase in population size by clutch and ratio of females
    for popu in 1:pops
        popsize1[:,:,popu] = [trunc(Int, popsize1[:,:,popu][1] * 2 *
                                        (1 - sex_bias[:,:,popu][1]) * clutch)]
    end


    #limit popsize to some maximum
    popsize1[popsize1 .> max_pop] .= max_pop

	#script crashes if pop = 0
	popsize1[popsize1 .<= 0] .= 1
    # #record the new popsize for plotting
    # total_pop = vcat(total_pop, popsize)


###################################################################
#Drift
###################################################################
    #drift
    for popu in 1:pops
        freq1[:,:,popu] = new_freq.(freq1[:,:,popu], popu)
		if sum(freq1[:,:,popu]) == 0
			freq1[1,:,popu] .= 1
		end
    end


    #sum freq to 1
    freq1 = mapslices(x -> x / sum(x), freq1, dims = [1,2])



###################################################################
#adjust arrays so that they dont include freq == 0
###################################################################

    # #find all haplotypes with freq != 0 which are retained
    # keep = freq1 .!= 0
    #
    # #only keep freqs with keep
    # freq1 = freq1[keep]
    #
    # #only keep the haplotypes with keep
    # pop1 = pop1[:,keep]
    #
    # # #only keep haplo fitness with keep
    # # haplo_fitness = haplo_fitness[:,keep]
    #
    # #only keep traits with keep
    # traits = traits[:,keep]
    #
    # #only keep drivers with keep
    # drivers = drivers[:,keep]
    #
    # #remove haplotypes that died in all
    # all = all[keep,:]


return(freq1, pop1, popsize1, drivers1, traits1, m_conver, sexes)



end


###################################################################
###################################################################
###################################################################
#running model as a single function
###################################################################
###################################################################
###################################################################


#
