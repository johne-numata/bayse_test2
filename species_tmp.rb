#! ruby -Ks
require 'csv'
require './thinkbayes'
require "numo/gnuplot"
require 'numo/narray'

class Subject
    # Represents a subject from the belly button study.
    attr_reader :code, :num_reads, :num_species, :total_reads, :total_species
	def initialize(code)
        # code: string ID
        # species: sequence of (int count, string species) pairs
        @code = code
        @species = []
        @suite = nil
        @num_reads = 0
        @num_species = 0
        @total_reads = 0
        @total_species = 0
        @prev_unseen = nil
        @pmf_n = nil
        @pmf_q = nil
        @pmf_l = nil
	end

    def add(species, count)
		# Add a species-count pair.
		# It is up to the caller to ensure that species names are unique.
		#   species: string species/genus name
        #   count: int number of individuals
        @species << [count, species]
	end

    def done(reverse = false, clean_param = 0)
        # Called when we are done adding species counts.
		#   reverse: which order to sort in
        clean(clean_param)  if clean_param != 0
        @species.reverse!	if reverse
        counts = get_counts()
        @num_species = counts.size
        @num_reads = counts.sum
    end

    def clean(clean_param = 50)
        # Identifies and removes bogus data.
		#   clean_param: parameter that controls the number of legit species
        def prob_bogus(k, r, clean_param)
            # Compute the probability that a species is bogus.
            q = clean_param / r
            p = (1 - q) ** k
            return p
		end
        print @code, clean_param

        counts = get_counts()
        r = 1.0 * counts.sum

        species_seq = []
        @species.sort.each{|k, species|
            next if rand() < prob_bogus(k, r, clean_param)
            species_seq << [k, species]
        }
        @species = species_seq
	end

    def get_m()
        # Gets number of observed species.
        return @species.size
    end

    def get_counts()
        # Gets the list of species counts
		# Should be in increasing order, if Sort() has been invoked.
        return @species.collect{|count, sp| count }
	end

    def make_cdf()
        # Makes a CDF of total prevalence vs rank.
        counts = get_counts()
        counts.reverse!
        new_counts = []
        counts.each_index{|i| new_counts << [i, counts[i]] }
        cdf = make_cdf_from_items(new_counts)
        return cdf
	end

    def get_names()
        # Gets the names of the seen species.
        return @species.collect{|c, name| name }
    end

    def get_species(index)
        # Gets the count and name of the indicated species.
		#   Returns: count-species pair
        return @species[index]
	end

    def process(low = nil, high: 500, conc: 1, iters: 100)
        # Computes the posterior distribution of n and the prevalences.
		#    low: minimum number of species
        #    high: maximum number of species
        #    conc: concentration parameter
        #    iters: number of iterations to use in the estimator
        counts = self.get_counts()
        m = counts.size
        low = [m, 2].max  if low == nil
        ns = low.step(high+1).to_a

        @suite = Species2.new(ns, conc, iters)
        @suite.update(counts)
	end

    def make_prediction(num_sims = 100)
        # Make predictions for the given subject.
		#   Precondition: Process has run
		#   num_sims: how many simulations to run for predictions
		#   pmf_l: predictive distribution of additional species
        add_reads = @total_reads - @num_reads
        curves = run_simulations(num_sims, add_reads)
        @pmf_l = make_predictive(curves)
	end

    def make_quick_prediction(num_sims = 100)
        # Make predictions for the given subject.
		#   Precondition: Process has run
		#   num_sims: how many simulations to run for predictions
		#   pmf_l: predictive distribution of additional species
        add_reads = @total_reads - @num_reads
        pmf = PmfSuite.new()
        temp, seen = get_seen_species()

        num_sims.times{
            temp, observations = generate_observations(add_reads)
            all_seen = seen.union(observations)
            l = all_seen.size - seen.size
            pmf.incr(l)
		}
        pmf.normalize()
        @pmf_l = pmf
	end

    def dist_l()
        # Returns the distribution of additional species, l.
        return @pmf_l
    end

    def make_figures()
        # Makes figures showing distribution of n and the prevalences.
        plot_dist_n()
        plot_prevalences()
	end

    def plot_dist_n()
        # Plots distribution of n.
        pmf = @suite.dist_n()
        print '90% CI for N: ', pmf.credible_interval(90), "\n"
        pmf.make_plot('Number of species', "species-ndist-#{@code}")
	end

    def plot_prevalences(num = 5)
        # Plots dist of prevalence for several species.
		#   num: how many species (starting with the highest prevalence)
        1.step(num + 1){|rank| plot_prevalence(rank) }
#                    axis=[0, 0.3, 0, 1],
     end

    def plot_prevalence(rank = 1, cdf_flag = true)
        # Plots dist of prevalence for one species.
		#   rank: rank order of the species to plot.
        #   cdf_flag: whether to plot the CDF
        # convert rank to index
        index = get_m() - rank

        temp, mix = @suite.dist_of_prevalence(index)
        count, temp = get_species(index)
#        mix.name = '%d (%d)' % (rank, count)

        print "90%% CI for prevalence of species #{rank}:" 
        print mix.credible_interval(90), "\n"

        if cdf_flag
            cdf = mix.make_cdf
            cdf.make_plot(nil ,"#{rank}")
        else
            mix.make_plot(nil ,"#{rank}")
		end
    end

    def get_seen_species()
        # Makes a set of the names of seen species.
		#   Returns: number of species, set of string species names
        names = get_names()
        m = names.size
        seen = species_generator(names, m).uniq
        return m, seen
	end

    def generate_observations(num_reads)
        # Generates a series of random observations.
		#    num_reads: number of reads to generate
		#    Returns: number of species, sequence of string species names
        n, prevalences = @suite.sample_posterior()

        names = get_names()
        name_iter = species_generator(names, n)

        items = name_iter.zip(prevalences)

        cdf = make_cdf_from_items(items)
        observations = cdf.sample(num_reads)

        #for ob in observations:
        #    print ob

        return n, observations
	end

    def run_simulation(num_reads, frac_flag = False, jitter = 0.01)
        # Simulates additional observations and returns a rarefaction curve.
		# k is the number of additional observations
        # num_new is the number of new species seen
		#    num_reads: how many new reads to simulate
        #    frac_flag: whether to convert to fraction of species seen
        #    jitter: size of jitter added if frac_flag is true
		#    Returns: list of (k, num_new) pairs
        m, seen = get_seen_species()
        n, observations = generate_observations(num_reads)

        curve = []
		observations.each_index{|i|
            seen.add(observation[i])

            if frac_flag
                frac_seen = len(seen) / float(n)
                frac_seen += rand(-jitter..jitter)
                curve << [i+1, frac_seen]
            else
                num_new = seen.size - m
                curve << [i+1, num_new]
            end
		}
        return curve
	end

    def run_simulations(num_sims, num_reads, frac_flag = false)
        # Runs simulations and returns a list of curves.
		# Each curve is a sequence of (k, num_new) pairs.
		#    num_sims: how many simulations to run
        #    num_reads: how many samples to generate in each simulation
        #    frac_flag: whether to convert num_new to fraction of total
        curves = []
        num_sims.times{ curves << run_simulation(num_reads, frac_flag) } 
        return curves
	end
end

def make_conditionals(curves, ks)
    # Makes Cdfs of the distribution of num_new conditioned on k.
	#   curves: list of (k, num_new) curves 
    #   ks: list of values of k
	#   Returns: list of Cdfs
    joint = make_joint_predictive(curves)
    cdfs = []
    ks.each{|k|
        pmf = joint.conditional(1, 0, k)
#        pmf.name = 'k=%d' % k
        cdf = pmf.make_cdf()
        cdfs << cdf
        print "90%% credible interval for #{k} \n"
        print cdf.credible_interval(90)
    }
    return cdfs
end

def make_joint_predictive(curves)
    # Makes a joint distribution of k and num_new.
	#   curves: list of (k, num_new) curves 
	#   Returns: joint Pmf of (k, num_new)
    joint = Joint.new()
    curves.each{|curve|  curve.each{|k, num_new|  joint.incr([k, num_new]) }}
    joint.normalize()
    return joint
end

def make_frac_cdfs(curves, ks)
    # Makes Cdfs of the fraction of species seen.
	#    curves: list of (k, num_new) curves 
	#    Returns: list of Cdfs
    d = {}
    curves.each{|curve|  curve.each{|k, frac| d[k] << frac  if ks.include(k) }}
    cdfs = {}
    d.items.ach{|k, fracs|
        cdf = make_cdf_from_list(fracs)
        cdfs[k] = cdf
	}
    return cdfs
end

## 以下微妙
def species_generator(names, num)
    # Generates a series of names, starting with the given names.
	# Additional names are 'unseen' plus a serial number.
	#    names: list of strings
    #    num: total number of species names to generate
	#    Returns: string iterator
    i = 0
    names.each{|name|
        yield name
        i += 1
	}
    while i < num do
        yield 'unseen-%d' % i
        i += 1
    end
end

def read_rarefacted_data(filename = 'journal.pone.0047712.s001.csv', clean_param = 0)
    # Reads a data file and returns a list of Subjects.
	# Data from http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0047712#s4
	#    filename: string filename to read
    #    clean_param: parameter passed to Clean
	#    Returns: map from code to Subject
    reader = CSV.read(filename)
    reader.shift
    
    subject = Subject.new('')
    subject_map = {}

    i = 0
    species = {}
    reader.each{|t|
        code = t[0]
        if code != subject.code
            subject = Subject.new(code)
            subject_map[code] = subject
		end
        # append a number to the species names so they're unique
        if !species.include?(t[1])
        	species[t[1]] = "#{t[1]}-#{i}"
        	i += 1
		end
        count = t[2].to_i
        subject.add(species[t[1]], count)
	}
    subject_map.each{|code, subject|   subject.done(clean_param) }
    return subject_map
end

def offset_curve(curve, i, n, dx = 0.3, dy = 0.3)
    # Adds random noise to the pairs in a curve.
	#   i is the index of the curve
    #   n is the number of curves
	# dx and dy control the amplitude of the noise in each dimension.
    xoff = -dx + 2 * dx * i / (n-1)
    yoff = -dy + 2 * dy * i / (n-1)
    curve = curve.collect{|x, y| [x + xoff, y + yoff] }
    return curve
end

def plot_curves(curves, root = 'species-rare')
    # Plots a set of curves.
	#   curves is a list of curves; each curve is a list of (x, y) pairs.
    n = curves.size
    curves.each_index{|i|
        curve = offset_curve(curves[i], i, n)
#未修正
        xs, ys = zip(*curve)
#		Numo.gnuplot do
#			set title:root
#			set xlabel:'# samples'
#  			set ylabel:'# species'
#			plot xs, ys, w:'lines'
#			sleep (2)
#		end
        thinkplot.Plot(xs, ys, color=color, alpha=0.3, linewidth=0.5)
	}
end

def plot_conditionals(cdfs, root = 'species-cond')
    # Plots cdfs of num_new conditioned on k.
	#   cdfs: list of Cdf
    #   root: string filename root

    cdfs.each{|cdf|  cdf.make_plot('# new species', root) }
end

def plot_frac_cdfs(cdfs, root = 'species-frac')
    # Plots CDFs of the fraction of species seen.
	#    cdfs: map from k to CDF of fraction of species seen after k samples

    cdfs.items.each{|k, cdf|
        xs, ys = cdf.render()
        ys = ys.collect{|y| 1 - y }
		Numo.gnuplot do
			set title:root
			set xlabel:'Fraction of species seen'
  			set ylabel:'Probability'
			plot xs, ys, w:'lines'
			sleep (2)
		end
	}
end

class Species < PmfSuite
    # Represents hypotheses about the number of species.
    def initialize(ns, conc = 1, iters = 1000)
        hypos = ns.collect{|n| Dirichlet.new(n, conc) }
        super(hypos)
        @iters = iters
    end

    def update(data)
        # Updates the suite based on the data.
		#     data: list of observed frequencies
        super(data)
        # update the next level of the hierarchy
        values.each{|hypo| hypo.update(data) }
	end

    def likelihood(data, hypo)
        # Computes the likelihood of the data under this hypothesis.
		#    hypo: Dirichlet object
        #    data: list of observed frequencies
        dirichlet = hypo

        # draw sample Likelihoods from the hypothetical Dirichlet dist
        # and add them up
        like = 0
        @iters.times{ like += dirichlet.likelihood(data) }

        # correct for the number of ways the observed species
        # might have been chosen from all species
        m = data.size
        like *= binomial_coef(dirichlet.n, m)
        return like
	end

    def dist_n()
        # Computes the distribution of n.
        pmf = PmfSuite.new()
        items.each{|hypo, prob| pmf.set(hypo.n, prob) }
        return pmf
     end
end 
 
class Species2
    # Represents hypotheses about the number of species.
	# Combines two layers of the hierarchy into one object.
	# ns and probs represent the distribution of N
	# params represents the parameters of the Dirichlet distributions    
    def initialize(ns, conc = 1, iters = 1000)
        @ns = ns
        @conc = conc
        @probs = Numo::DFloat.ones(ns.size)
        @params = Numo::DFloat.ones(@ns[-1]) * conc
        @iters = iters
        @num_reads = 0
        @m = 0
	end

    def update(data)
        # Updates the distribution based on data.
		#    data: numpy array of counts
        @num_reads += data.sum

        like = Numo::DFloat.zeros(@ns.size)
        @iters.times{ like += sample_likelihood(data) }

        @probs *= like
        @probs /= @probs.sum

        @m = data.size
        #self.params[:self.m] += data * self.conc   # python元ファイルでのコメント
        @params[0..@m-1] += data
	end

    def sample_likelihood(data)
        # Computes the likelihood of the data for all values of n.
		# Draws one sample from the distribution of prevalences.
		#   data: sequence of observed counts
		#   Returns: numpy array of m likelihoods
        gammas = @params.to_a.collect{|r| Random::DEFAULT.gamma(r) }

        m = data.size
        row = Numo::DFloat.new(m)
        row.store(gammas[0..m-1])
        col = gammas.size.times.map{|i| gammas[0..i].inject(:+)}

        log_likes = []
        @ns.each{|n|
            ps = row / col[n - 1]
            terms = Numo::NMath.log(ps) * data
            log_like = terms.sum()
            log_likes << log_like
		}
        max_l = log_likes.max
        log_likes.collect{|a| a -= max_l } 
        likes = log_likes.collect{|a| Numo::NMath.exp(a) }
#        log_likes -= log_likes.max
#        likes = Math.exp(log_likes)

#        coefs = @ns.collect{|n| binomial_coef(n, m) }
#        likes *= coefs
#        return likes
        coefs = @ns.collect{|n| binomial_coef(n, m) }
        return likes.collect{|a| a *= coefs }
	end

    def dist_n()
        # Computes the distribution of n.
		#    Returns: new Pmf object
        pmf = make_pmf_from_items(@ns.zip(@probs))
        return pmf
	end

    def random_n()
        # Returns a random value of n.
        return dist_n().random()
	end

    def marginal_beta(n, index)
        # Computes the conditional distribution of the indicated species.
        #   n: conditional number of species
        #   index: which species
		#   Returns: Beta object representing a distribution of prevalence.
        alpha0 = @params[0..n-1].sum
        alpha = @params[index]
        return Beta.new(alpha, alpha0 - alpha)
	end

    def dist_of_prevalence(index)
        # Computes the distribution of prevalence for the indicated species.
		#    index: which species
		#    Returns: (metapmf, mix) where metapmf is a MetaPmf and mix is a Pmf
        metapmf = PmfSuite.new()

        @ns.zip(@probs.to_a).each{|n, prob|
            beta = marginal_beta(n, index)
            pmf = beta.make_pmf()
            metapmf.set(pmf, prob)
		}
        mix = make_mixture(metapmf)
        return metapmf, mix
    end

    def sample_posterior()
        # Draws random n and prevalences.
		#   Returns: (n, prevalences)
        n = random_n()
        prevalences = sample_prevalences(n)

        #print 'Peeking at n_cheat'
        #n = n_cheat

        return n, prevalences
	end

    def sample_prevalences(n)
        # Draws a sample of prevalences given n.
		#    n: the number of species assumed in the conditional
		#    Returns: numpy array of n prevalences
        return [1.0]  if n == 1

        q_desired = random_q(n)
        q_desired = Math.max(q_desired, 1e-6)

        params = unbias(n, @m, q_desired)

        gammas = rand_gamma(params)
        gammas /= gammas.sum()
        return gammas
    end

    def unbias(n, m, q_desired)
        # Adjusts the parameters to achieve desired prev_unseen (q).
		#    n: number of species
        #    m: seen species
        #    q_desired: prevalence of unseen species
        params = @params[0..n-1].clone

        return params  if n == m
        
        x = params[0..m-1].sum
        y = params[m - 1 .. -1].sum
        a = x + y
        #print x, y, a, x/a, y/a

        g = q_desired * a / y
        f = (a - g * y) / x
        params.each_index{|i|
        	if i < m
        		params[i] *= f
        	else
        		params[i] *= g
			end
		}
        return params
	end
end

def simple_dirichlet_example()
    # Makes a plot showing posterior distributions for three species.
	# This is the case where we know there are exactly three species.
    names = ['lions',  'tigers', 'bears']
    data = [3, 2, 1]

    dirichlet = Dirichlet.new(3)
    0.step(2){|i|
        beta = dirichlet.marginal_beta(i)
        beta.name = names[i]
        print 'mean = ', names[i], " ", beta.mean(), "\n"
	}
    dirichlet.update(data)
    t = 0.step(2).collect{|i|
        beta = dirichlet.marginal_beta(i)
        beta.name = names[i]
        print 'mean = ', names[i], " ", beta.mean(), "\n"
        pmf = beta.make_pmf
        pmf.render
	}
	make_plot_2d(t, x_label: 'Prevalence', y_label: "PMF", title: 'species1', grid: "on", )
end

def hierarchical_example()
    # Shows the posterior distribution of n for lions, tigers and bears.
    ns = 3.step(30).to_a
    suite = Species.new(ns, 1.0)

    data = [3, 2, 1]
    suite.update(data)

    pmf = suite.dist_n()
    pmf.make_plot(x_label:'Number of species', y_label:'PMF', title:'species')
end

def run_subject(code, conc = 1, high = 500)
    # Run the analysis for the subject with the given code.
	#   code: string code
    subjects = read_rarefacted_data()
    subject = subjects[code]

    subject.process(conc: conc, high: high, iters: 300)
#    subject.make_prediction()
#    subject.make_quick_prediction()

    print_summary(subject)
    actual_l = subject.total_species - subject.num_species
    cdf_l = subject.dist_l().make_cdf()
    print_prediction(cdf_l, actual_l)

    subject.make_figures()

    num_reads = 400
    curves = subject.run_simulations(100, num_reads)
    root = "species-rare-#{subject.code}"
    plot_curves(curves, root)

    num_reads = 800
    curves = subject.run_simulations(500, num_reads)
    ks = [100, 200, 400, 800]
    cdfs = make_conditionals(curves, ks)
    root = "species-cond-#{subject.code}"
    plot_conditionals(cdfs, root)

    num_reads = 1000
    curves = subject.run_simulations(500, num_reads, frac_flag = true)
    ks = [10, 100, 200, 400, 600, 800, 1000]
    cdfs = make_frac_cdfs(curves, ks)
    root = "species-frac-#{subject.code}"
    Plot_frac_cdfs(cdfs, root)
end

def print_summary(subject)
    # Print a summary of a subject.
	#     subject: Subject
    print subject.code
    print " found #{subject.num_species} species in #{subject.num_reads} reads\n" 
    print " total #{subject.total_species} species in #{subject.total_reads} reads\n"
    cdf = subject.suite.dist_n().make_cdf()
    print "n\n"
    print_prediction(cdf, 'unknown')
end

def print_prediction(cdf, actual)
    # Print a summary of a prediction.
	#    cdf: predictive distribution
    #    actual: actual value
    median = cdf.percentile(50)
    low, high = cdf.credible_interval(75)
    printf "predicted %0.2f (%0.2f %0.2f)\n", median, low, high
    print "actual = ", actual, "\n"
end

####################################################
###   メインルーチン
####################################################
 run_subject('B1242', conc = 1, high = 100)

#  simple_dirichlet_example()

#  hierarchical_example()
