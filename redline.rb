#! ruby -Ks
require 'csv'
require './thinkbayes'
require './kde'
require "numo/gnuplot"
require "numo/narray"

# z: time between trains
# x: time since the last train
# y: time until the next train
# zb: distribution of z as seen by a random arrival
# longest hypothetical time between trains, in seconds

UPPER_BOUND = 1200

# observed gaps between trains, in seconds
# collected using code in redline_data.py, run daily 4-6pm
# for 5 days, Monday 6 May 2013 to Friday 10 May 2013
OBSERVED_GAP_TIMES = [
    428.0, 705.0, 407.0, 465.0, 433.0, 425.0, 204.0, 506.0, 143.0, 351.0, 
    450.0, 598.0, 464.0, 749.0, 341.0, 586.0, 754.0, 256.0, 378.0, 435.0, 
    176.0, 405.0, 360.0, 519.0, 648.0, 374.0, 483.0, 537.0, 578.0, 534.0, 
    577.0, 619.0, 538.0, 331.0, 186.0, 629.0, 193.0, 360.0, 660.0, 484.0, 
    512.0, 315.0, 457.0, 404.0, 740.0, 388.0, 357.0, 485.0, 567.0, 160.0, 
    428.0, 387.0, 901.0, 187.0, 622.0, 616.0, 585.0, 474.0, 442.0, 499.0, 
    437.0, 620.0, 351.0, 286.0, 373.0, 232.0, 393.0, 745.0, 636.0, 758.0,
]

def bias_pmf(pmf, invert = false)
    # Returns the Pmf with oversampling proportional to value.
    # If invert=True, computes in inverse operation; for example,
    # unbiasing a sample collected from students.
    new_pmf = pmf.copy()
    for x in pmf.values()
        if invert
            new_pmf.mult(x, 1.0/x)
        else
            new_pmf.mult(x, x)
        end
    end
    new_pmf.normalize()
    return new_pmf
end

def unbias_pmf(pmf)
    # Returns the Pmf with oversampling proportional to 1/value.
    return bias_pmf(pmf, true)
end

def make_uniform_pmf(low, high)
    # Make a uniform Pmf.
	#   low: lowest value (inclusive)
    #   high: highest value (inclusive)
    pmf = PmfSuite.new()
    for x in make_range(low = low, high = high)
        pmf.set(x, 1)
    end
    pmf.normalize()
    return pmf    
end

def make_range(low = 10, high = nil, skip = 10)
    # Makes a range representing possible gap times in seconds.
    high = UPPER_BOUND   if high.nil?
    return low.step(high + skip, skip).to_a
end

class WaitTimeCalculator
    # Encapsulates the forward inference process.
	# 
    # Given the actual distribution of gap times (z),
    # computes the distribution of gaps as seen by
    # a random passenger (zb), which yields the distribution
    # of wait times (y) and the distribution of elapsed times (x).
	attr_reader :pmf_x, :pmf_zb
	def initialize(pmf, inverse = false)
	 	# pmf: Pmf of either z or zb
        # inverse: boolean, true if pmf is zb, false if pmf is z
        if inverse
            @pmf_zb = pmf
            @pmf_z = unbias_pmf(pmf)
        else
            @pmf_z = pmf
            @pmf_zb = bias_pmf(pmf)
        end
        # distribution of wait time
        @pmf_y = pmf_of_wait_time(@pmf_zb)

        # the distribution of elapsed time is the same as the
        # distribution of wait time
        @pmf_x = @pmf_y
    end
	
    def generate_sample_wait_times(n)
        # Generates a random sample of wait times.
        cdf_y = make_cdf_from_pmf(@pmf_y)
        sample = cdf_y.sample(n)
        return sample
    end
	
    def generate_sample_gaps(n)
		# Generates a random sample of gaps seen by passengers.
        cdf_zb = make_cdf_from_pmf(@pmf_zb)
        sample = cdf_zb.sample(n)
        return sample
    end
	
    def generate_sample_passengers(lam, n)
        # Generates a sample wait time and number of arrivals.
		#   lam: arrival rate in passengers per second
        #   n: number of samples
		#   Returns: list of (k1, y, k2) tuples
        #     k1: passengers there on arrival
        #     y: wait time
        #     k2: passengers arrived while waiting
        zs = generate_sample_gaps(n)
        xs, ys = split_gaps(zs)

        res = []
        xs.zip(ys).each{| x, y |
            k1 = Random::DEFAULT.poisson(x * lam)
            k2 = Random::DEFAULT.poisson(y * lam)
            res << [k1, y, k2]
		}
        return res
    end

    def plot_pmfs(root = 'redline0')
     	# Plots the computed Pmfs.
	 	# root: string
        pmfs = scale_dists([@pmf_z, @pmf_zb], 1.0/60)
	    z = pmfs[0].items.transpose;  zb = pmfs[1].items.transpose
		t = [[z[0], z[1], w:'lines', t:"z"], [zb[0], zb[1], w:'lines', t:"zb"]]
		make_plot_2d(t, x_label: "Time (min)", y_label: "CDF", title: root, grid: "on", )
    end

    def make_plot(root='redline2')
        # Plots the computed CDFs.
	 	# root: string
        print 'Mean z', @pmf_z.mean() / 60
        print 'Mean zb', @pmf_zb.mean() / 60
        print 'Mean y', @pmf_y.mean() / 60
        print "\n"

        cdf_z = @pmf_z.make_cdf()
        cdf_zb = @pmf_zb.make_cdf()
        cdf_y = @pmf_y.make_cdf()

        cdfs = scale_dists([cdf_z, cdf_zb, cdf_y], 1.0/60)
	    z = cdfs[0].items.transpose
	    zb = cdfs[1].items.transpose
	    y = cdfs[2].items.transpose
		t = [[z[0], z[1], w:'lines', t:"Z"], [zb[0], zb[1], w:'lines', t:"zb"],
				 [y[0], y[1], w:'lines', t:"y"]]
		make_plot_2d(t, x_label: "Time (min)", y_label: "CDF", title: root, grid: "on", )
    end
end

def split_gaps(zs)
    # Splits zs into xs and ys.
	# zs: sequence of gaps
	# Returns: tuple of sequences (xs, ys)

    xs = zs.collect{|z| rand(z)}
    ys = zs.zip(xs).collect{|z, x| z - x}
    return xs, ys
end

def pmf_of_wait_time(pmf_zb)
    # Distribution of wait time.
	# pmf_zb: dist of gap time as seen by a random observer
	# Returns: dist of wait time (also dist of elapsed time)
    metapmf = PmfSuite.new()
    pmf_zb.items.each{| gap, prob |
        uniform = make_uniform_pmf(0, gap)
        metapmf.set(uniform, prob)
    }
    pmf_y = make_mixture(metapmf)
    return pmf_y
end

def scale_dists(dists, factor)
    # Scales each of the distributions in a sequence.
    # dists: sequence of Pmf or Cdf
    # factor: float scale factor
    return dists.collect{|dist| dist.scale(factor)}
end

class ElapsedTimeEstimator
    # Uses the number of passengers to estimate time since last train."""
	attr_reader :pmf_y, :post_x

    def initialize(wtc, lam, num_passengers)
		# pmf_x: expected distribution of elapsed time
		# lam: arrival rate in passengers per second
        # num_passengers: # passengers seen on the platform

        # prior for elapsed time
        @prior_x = Elapsed.new(wtc.pmf_x)

        # posterior of elapsed time (based on number of passengers)
        @post_x = Marshal.load(Marshal.dump(@prior_x))
        @post_x.update([lam, num_passengers])

        # predictive distribution of wait time
        @pmf_y = predict_wait_time(wtc.pmf_zb, @post_x)
    end

    def make_plot()
        # observed gaps
        cdf_prior_x = @prior_x.make_cdf()
        cdf_post_x = @post_x.make_cdf()
        cdf_y = @pmf_y.make_cdf()

        cdfs = scale_dists([cdf_prior_x, cdf_post_x, cdf_y], 1.0/60)
	   	z = cdfs[0].items.transpose; zb = cdfs[1].items.transpose; y = cdfs[2].items.transpose
		t = [[z[0], z[1], w:'lines', t:"prior_x"], [zb[0], zb[1], w:'lines', t:"poterior_x"],
				 [y[0], y[1], w:'lines', t:"y"]]
		make_plot_2d(t, x_label: "Time (min)", y_label: "CDF", title: "åoâﬂéûä‘ó\ë™", grid: "on", )
    end
end

class ArrivalRate < PmfSuite
    # Represents the distribution of arrival rates (lambda)."""

    def likelihood(data, hypo)
 		# Computes the likelihood of the data under the hypothesis.
		# Evaluates the Poisson PMF for lambda and k.
		# hypo: arrival rate in passengers per second
        # data: tuple of elapsed_time and number of passengers

        lam = hypo
        x, k = data
        like = eval_poisson_pmf(k, lam * x)
        return like
	end
end

class ArrivalRateEstimator
    # Estimates arrival rate based on passengers that arrive while waiting.
    attr_reader  :post_lam

    def initialize(passenger_data)
    	# passenger_data: sequence of (k1, y, k2) pairs

        low, high = 0, 5
        n = 51
        hypos = low.step(high, (high.to_f - low.to_f)/(n - 1)).to_a.collect{|a| a / 60}

        @prior_lam = ArrivalRate.new(hypos)
        @prior_lam.remove(0.0)

        @post_lam = @prior_lam.copy()

        for _k1, y, k2 in passenger_data
            @post_lam.update([y, k2])
		end
        print 'Mean posterior lambda', @post_lam.mean(), "\n"
	end
	
    def make_plot(root = "redline1")
        # Plot the prior and posterior CDF of passengers arrival rate.
        # convert units to passengers per minute
        prior = @prior_lam.make_cdf().scale(60).items.transpose
        post = @post_lam.make_cdf().scale(60).items.transpose
		make_plot_2d([[prior[0], prior[1], w:'lines', t: "prior"],[post[0], post[1], w:'lines', t:"posterior"]],
						 x_label: "Arrival rate (passengers / min)", y_label: "CDF", title: root, grid: "on", )
	end
end

class Elapsed < PmfSuite
    # Represents the distribution of elapsed time (x)."""

    def likelihood(data, hypo)
        # Evaluates the Poisson PMF for lambda and k.
		# hypo: elapsed time since the last train
        # data: tuple of arrival rate and number of passengers

        x = hypo
        lam, k = data
        like = eval_poisson_pmf(k, lam * x)
        return like
    end
end

def predict_wait_time(pmf_zb, pmf_x)
    # Computes the distribution of wait times.
	# Enumerate all pairs of zb from pmf_zb and x from pmf_x,
    # and accumulate the distribution of y = z - x.
	# pmf_zb: distribution of gaps seen by random observer
    # pmf_x: distribution of elapsed time

    pmf_y = pmf_zb - pmf_x
    remove_negatives(pmf_y)
    return pmf_y
end

def remove_negatives(pmf)
    # Removes negative values from a PMF.
	# pmf: Pmf

    pmf.values().each{|val|  pmf.remove(val)  if val < 0}
    pmf.normalize()
end

class Gaps < PmfSuite
    # Represents the distribution of gap times,
    # as updated by an observed waiting time.
    def likelihood(data, hypo)
		# hypo: actual time between trains
        # data: observed wait time
        z = hypo
        y = data
        return 0  if y > z
        return 1.0 / z
	end
end

class WaitMixtureEstimator
    # Encapsulates the process of estimating wait time with uncertain lam.

    def initialize(wtc, are, num_passengers = 15)
		# wtc: WaitTimeCalculator
        # are: ArrivalTimeEstimator
        # num_passengers: number of passengers seen on the platform

        @metapmf = PmfSuite.new()

        are.post_lam.items().sort.each{|lam, prob|
            ete = ElapsedTimeEstimator.new(wtc, lam, num_passengers)
            @metapmf.set(ete.pmf_y, prob)
		}
        @mixture = make_mixture(@metapmf)

        lam = are.post_lam.mean()
        ete = ElapsedTimeEstimator.new(wtc, lam, num_passengers)
        @point = ete.pmf_y
	end

    def make_plot(root = "redline4")
        # plot the MetaPmf
		t = []
        @metapmf.items.each{| pmf, prob |
        	cdf = pmf.make_cdf().scale(1.0/60)
	        pro = cdf.items.transpose
			t << [pro[0], pro[1], w:'lines']
		}
		make_plot_2d(t, x_label: "ë“Çøéûä‘", y_label: "CDF", title: root, grid: "on", key: nil)

        # plot the mixture and the distribution based on a point estimate
        po = @point.make_cdf.scale(1.0/60).items.transpose
        mix = @mixture.make_cdf.scale(1.0/60).items.transpose
		make_plot_2d([[po[0], po[1], w:'lines', t: "point"],[mix[0], mix[1], w:'lines', t: "mixture"]],
						 x_label: "ë“Çøéûä‘", y_label: "CDF", title: root, grid: "on", )
	end
end


def generate_sample_data(gap_times, lam=0.0333, n=10)
    # Generates passenger data based on actual gap times.
	# gap_times: sequence of float
    # lam: arrival rate in passengers per second
    # n: number of simulated observations

    xs = make_range(low=10)
    pdf_z = estimated_pdf(gap_times)
    pmf_z = pdf_z.make_pmf(xs)

    wtc = wait_time_calculator(pmf_z, inverse = False)
    passenger_data = wtc.generate_sample_passengers(lam, n)
    return wtc, passenger_data
end

def run_simple_process(gap_times, lam = 0.0333, num_passengers = 15, plot = true)
    # Runs the basic analysis and generates figures.
	# gap_times: sequence of float
    # lam: arrival rate in passengers per second
    # num_passengers: int number of passengers on the platform
    # plot: boolean, whether to generate plots
	# Returns: WaitTimeCalculator, ElapsedTimeEstimator

    cdf_z = make_cdf_from_list(gap_times).scale(1.0/60)
    print 'CI z = ', cdf_z.credible_interval(90), "\n"

    xs = make_range(low = 10)
    pdf_z = EstimatedPdf.new(gap_times)
    pmf_z = pdf_z.make_pmf(xs)

    wtc = WaitTimeCalculator.new(pmf_z, inverse = false)    
    if plot
        wtc.plot_pmfs()
        wtc.make_plot()
    end

    ete = ElapsedTimeEstimator.new(wtc, lam, num_passengers)
    if plot
        ete.make_plot()
	end
	
    return wtc, ete
end

def run_mix_process(gap_times, lam = 0.0333, num_passengers = 15, plot = true)
    # Runs the analysis for unknown lambda.
    #    gap_times: sequence of float
    #    lam: arrival rate in passengers per second
    #    num_passengers: int number of passengers on the platform
    #    plot: boolean, whether to generate plots
    #    Returns: WaitMixtureEstimator

    wtc, ete = run_simple_process(gap_times, lam, num_passengers)
    passenger_data = wtc.generate_sample_passengers(lam, n = 5)

    total_y = 0
    total_k2 = 0
    for k1, y, k2 in passenger_data
        print k1," ", y/60, " ", k2, "\n"
        total_y += y/60
        total_k2 += k2
    end
    print total_k2, " ", total_y, "\n"
    print "Average arrival rate  ", total_k2 / total_y, "\n"

    are = ArrivalRateEstimator.new(passenger_data)
    are.make_plot()   if plot

    wme = WaitMixtureEstimator.new(wtc, are, num_passengers)
    wme.make_plot() if plot

    return wme
end

def run_loop(gap_times, nums, lam = 0.0333)
    # Runs the basic analysis for a range of num_passengers.
	#   gap_times: sequence of float
    #   nums: sequence of values for num_passengers
    #   lam: arrival rate in passengers per second
	#   Returns: WaitMixtureEstimator
	
    # resample gap_times
    n = 220
    cdf_z = make_cdf_from_list(gap_times)
    sample_z = cdf_z.sample(n)
    pmf_z = make_suite_from_list(sample_z)

    # compute the biased pmf and add some long delays
    cdf_zp = bias_pmf(pmf_z).make_cdf()
    sample_zb = cdf_zp.sample(n) + [1800, 2400, 3000]

    # smooth the distribution of zb
    pdf_zb = EstimatedPdf.new(sample_zb)
    xs = make_range(low = 60, nil, 5)
    pmf_zb = pdf_zb.make_pmf(xs)

    # unbias the distribution of zb and make wtc
    pmf_z = unbias_pmf(pmf_zb)
    wtc = WaitTimeCalculator.new(pmf_z)

    probs = []; t = []
    for num_passengers in nums
        ete = ElapsedTimeEstimator.new(wtc, lam, num_passengers)

        # compute the posterior prob of waiting more than 15 minutes
        cdf_y = ete.pmf_y.make_cdf()
        probs << 1 - cdf_y.prob(900)
        pro = ete.pmf_y.make_cdf.items.transpose
       	t << [pro[0], pro[1], w:'lines', t: num_passengers.to_s]
    end
   	make_plot_2d(t, x_label: "ë“Çøéûä‘", y_label: "PMF", title: "redline5", grid: "on", )
   	make_plot_2d([nums, probs, w:'lines'], x_label: "Num passengers", y_label: "P(y > 15 min)", title: "redline5",
   						 logscale: "y",  grid: "on", )
end

#########################################################################
### ÉÅÉCÉìÉãÅ[É`Éì
#########################################################################
run_loop(OBSERVED_GAP_TIMES, nums = [0, 5, 10, 15, 20, 25, 30, 35])
run_mix_process(OBSERVED_GAP_TIMES)

#########################################################################
###   à»â∫ìdé‘ä‘äuÇÃêÑë™Ç…ä÷Ç∑ÇÈãcò_Ç…Ç¬Ç¢ÇƒÇÃñÕîÕâìö?
#########################################################################
class GapDirichlet  <  Dirichlet
    # Represents the distribution of prevalences for each gap time.
    def initialize(xs)
        # xs: sequence of possible gap times
        n = xs.size
        super(n)
        @xs = xs
        @mean_zbs = []
	end

    def pmf_mean_zb()
        # Makes the Pmf of mean zb.
		# Values stored in mean_zbs.

        return make_suite_from_list(@mean_zbs)
	end

    def preload(data)
        # Adds pseudocounts to the parameters.
		# data: sequence of pseudocounts

        public_method(:update).super_method.call (data)
	end

    def update(data)
        # Computes the likelihood of the data.
		# data: wait time observed by random arrival (y)
		# Returns: float probability

        k, y = data

        print k, y
        prior = predictive_pmf(@xs)
        gaps = Gaps.new(prior)
        gaps.update(y)
        probs = gaps.probs(@xs)

        @params.each_index{|i| @params[i] += probs[i]}
	end
end

class GapDirichlet2  <  GapDirichlet
    # Represents the distribution of prevalences for each
    # gap time.

    def update(data)
        # Computes the likelihood of the data.
		# data: wait time observed by random arrival (y)
		# Returns: float probability

        k, y = data

        # get the current best guess for pmf_z
        pmf_zb = predictive_pmf(@xs)

        # use it to compute prior pmf_x, pmf_y, pmf_z
        wtc = WaitTimeCalculator.new(pmf_zb, inverse = true)

        # use the observed passengers to estimate posterior pmf_x
        elapsed = ElapsedTimeEstimator.new(wtc,
                                       lam=0.0333,
                                       num_passengers=k)

        # use posterior_x and observed y to estimate observed z
        obs_zb = elapsed.post_x + floor(y)
        probs = obs_zb.probs(@xs)

        mean_zb = obs_zb.mean()
        @mean_zbs << mean_zb
        print k, y, mean_zb

        # use observed z to update beliefs about pmf_z
        @params.each_index{|i| @params[i] += probs[i] }
	end
end

class GapTimeEstimator
    # Infers gap times using passenger data.

    def initialize(xs, pcounts, passenger_data)
        @xs = xs
        @pcounts = pcounts
        @passenger_data = passenger_data

        @wait_times = passenger_data.collect{|_k1, y, _k2| y}
        @pmf_y = make_suite_from_list(@wait_times, name: "y")

        dirichlet = GapDirichlet2.new(@xs)
        dirichlet.params.collect{|a| a /= 1.0 }

        dirichlet.preload(@pcounts)
        dirichlet.params.collect{|a| a /= 20.0 }

        @prior_zb = dirichlet.predictive_pmf(@xs)
        
        for k1, y, _k2 in passenger_data
            dirichlet.update([k1, y])
		end
        @pmf_mean_zb = dirichlet.pmf_mean_zb()

        @post_zb = dirichlet.predictive_pmf(@xs); @post_zb.name = "zb"
        @post_z = unbias_pmf(@post_zb); @post_z.name = "z"
	end

    def plot_pmfs()
        print 'Mean y', @pmf_y.Mean()
        print 'Mean z', @post_z.Mean()
        print 'Mean zb', @post_zb.Mean()

		t = [@pmf_y.render, @post_z.render, @post_zb.render]
		make_plot_2d(t, x_label: "Time (min)", y_label: "PMF", title: "ìdé‘ä‘äuÇÃó\ë™", grid: "on", )
	end

    def make_plot()
		t = [@pmf_y.make_cdf.render, @proir_zb.make_cdf.render, @post_zb.make_cdf.render, @pmf_mean_zb.make_cdf.render]
		make_plot_2d(t, x_label: "Time (min)", y_label: "PMF", title: "ìdé‘ä‘äuÇÃó\ë™", grid: "on", )
	end
end

def floor(x, factor = 10)
    # Rounds down to the nearest multiple of factor.
	# When factor=10, all numbers from 10 to 19 get floored to 10.
    return (x / factor).to_i * factor
end

def test_gte()
    # Tests the GapTimeEstimator.
    xs = [60, 120, 240]
    gap_times = [60, 60, 60, 60, 60, 120, 120, 120, 240, 240]

    # distribution of gap time (z)
    pdf_z = EstimatedPdf.new(gap_times)
    pmf_z = pdf_z.make_pmf(xs)

    wtc = WaitTimeCalculator.new(pmf_z, inverse = false)

    lam = 0.0333
    n = 100
    passenger_data = wtc.generate_sample_passengers(lam, n)
    pcounts = [0, 0, 0]
    ite = GapTimeEstimator.new(xs, pcounts, passenger_data)
#	 thinkplot.Cdf(wtc.pmf_z.MakeCdf(name="actual z"))    
#    thinkplot.Cdf(wtc.pmf_zb.MakeCdf(name="actual zb"))
    ite.make_plot()
end

#GapTimeêÑë™ÇÕñ¢äÆê¨
#test_gte
