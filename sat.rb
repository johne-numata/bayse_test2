#! ruby -Ks
require './thinkbayes'
require "numo/gnuplot"
require 'csv'

def read_scale(filename = 'sat_scale.csv', col = 2)
    # Reads a CSV file of SAT scales (maps from raw score to standard score).
    #   col: which column to start with (0=Reading, 2=Math, 4=Writing)
    def parse_range(s)
        # Parse a range of values in the form 123-456
		#   s: string
        t = s.split('-').collect{|x| x.to_i }
        return 1.0 * t.sum / t.size
	end
	reader = CSV.read(filename)
    raws = []
    scores = []

    reader.each{|t|
        if t[col] && t[col + 1] && t[col].int_valid? && \
        		 (t[col + 1][0].int_valid? || t[col + 1][1].int_valid?)
            raws << t[col].to_i 
            scores << parse_range(t[col + 1])
        end
	}
    return Interpolator.new(raws.sort, scores.sort)  # 昇順ソートが必要
end

def read_ranks(filename = 'sat_ranks.csv')
	#   Returns: list of (score, freq) pairs
    reader = CSV.read(filename)
    res = []

    reader.each{|t|
        if t[0] && t[0].int_valid? && t[1].int_valid?
            score = t[0].to_i
            freq = t[1].to_i
            res << [score, freq]
        end
	}
    return res
end

def divide_values(pmf, denom)
    # Divides the values in a Pmf by denom.
	#    Returns a new Pmf
    new = PmfSuite.new()
    denom = denom.to_f
    pmf.items.each{|val, prob|
        x = val / denom
        new.set(x, prob)
    }
    return new
end

class Exam
    # Encapsulates information about an exam.
	#   Contains the distribution of scaled scores and an
    #   Interpolator that maps between scaled and raw scores.
    attr_reader :prior, :max_score
    def initialize
        @scale = read_scale()
        scores = read_ranks()
        score_pmf = make_suite_from_dict(scores.to_h)
        @raw = reverse_scale(score_pmf)
        @max_score = @raw.values().max
        @prior = divide_values(@raw, denom = @max_score)
        
        center = -0.05		# 実際の分布に合わせるための係数、他で試算した結果
        width = 1.8			# 実際の分布に合わせるための係数、他で試算した結果
        @difficulties = make_difficulties(center, width, @max_score)
	end

    def compare_scores(a_score, b_score, constructor)
        # Computes posteriors for two test scores and the likelihood ratio.
		#   a_score, b_score: scales SAT scores
        #   constructor: function that instantiates an Sat or Sat2 object
        if constructor == "Sat"
        	a_sat = Sat.new(self, a_score)
        	b_sat = Sat.new(self, b_score)
        elsif constructor == "Sat2"
        	a_sat = Sat2.new(self, a_score)
        	b_sat = Sat2.new(self, b_score)
        end

        a_sat.plot_posteriors(b_sat)

        if constructor == "Sat"
            plot_joint_dist(a_sat, b_sat)
		end
        top = TopLevel.new('AB')
        top.update([a_sat, b_sat])
        top.print_dict()

        ratio = top.prob('A') / top.prob('B')
        
        print 'Likelihood ratio = ', ratio, "\n"

        posterior = ratio / (ratio + 1)
        print 'Posterior = ', posterior, "\n"

        if constructor == "Sat2"
            compare_posterior_predictive(a_sat, b_sat)
        end
	end

    def make_raw_score_dist(efficacies)
        # Makes the distribution of raw scores for given difficulty.
		#   efficacies: Pmf of efficacy
        pmfs = PmfSuite.new()
        efficacies.items.each{|efficacy, prob|
            scores = pmf_correct(efficacy)
            pmfs.set(scores, prob)
		}
        mix = make_mixture(pmfs)
        return mix
	end

    def calibrate_difficulty
        # Make a plot showing the model distribution of raw scores.
        cdf = make_cdf_from_pmf(@raw)
		make_plot_2d(cdf.render, x_label: "raw score", y_label: "CDF", title: "sat_calibrate")

        efficacies = make_gaussian_pmf(0, 1.5, 3)
        pmf = make_raw_score_dist(efficacies)
        cdf = make_cdf_from_pmf(pmf)
		make_plot_2d(cdf.render, x_label: "raw score", y_label: "CDF", title: "sat_calibrate")
	end

    def pmf_correct(efficacy)
        # Returns the PMF of number of correct responses.
		#    efficacy: float
        pmf = toplevel_pmf_correct(efficacy, @difficulties)
        return pmf
	end

    def lookup(raw)
        # Looks up a raw score and returns a scaled score.
        return @scale.lookup(raw)
	end
        
    def reverse(score)
        # Looks up a scaled score and returns a raw score.
		# Since we ignore the penalty, negative scores round up to zero.
        raw = @scale.reverse(score)
        return raw > 0 ? raw : 0
	end
        
    def reverse_scale(pmf)
        # Applies the reverse scale to the values of a PMF.
		#   pmf: Pmf object
        #   scale: Interpolator object
		#   Returns: new Pmf
        new = PmfSuite.new()
        pmf.items.each{|val, prob|
            raw = reverse(val)
            new.incr(raw, prob)
        }
        return new
	end
end

class Sat < PmfSuite
    # Represents the distribution of p_correct for a test-taker.
    attr_reader :score
    def initialize(exam, score)
        @exam = exam
        @score = score

        # start with the prior distribution
        super (exam.prior)

        # update based on an exam score
        update(score)
    end

    def likelihood(data, hypo)
        # Computes the likelihood of a test score, given efficacy.
        p_correct = hypo
        score = data

        k = @exam.reverse(score)
        n = @exam.max_score
        like = eval_binomial_pmf(k, n, p_correct)
        return like
	end

    def plot_posteriors(other)
        # Plots posterior distributions of efficacy.
		#  other: Sat objects.
        cdf1 = make_cdf_from_pmf(self)
        cdf2 = make_cdf_from_pmf(other)
		d1 = cdf1.items.transpose
		d2 = cdf2.items.transpose
		Numo.gnuplot do
#			set axis: "0.7, 1.0, 0.0, 1.0"
       		set title: 'sat_posteriors_p_corr'
			set xlabel: 'p_correct'
			set ylabel: "cdf"
		    plot [d1[0], d1[1], w:'lines', t:'self'], [d2[0], d2[1], w:'lines', t:'other']
		    sleep (2)
    	end
 	end
end

class Sat2 < PmfSuite
    # Represents the distribution of efficacy for a test-taker.
    attr_reader :score
    def initialize(exam, score)
        @exam = exam
        @score = score
        # start with the Gaussian prior
        efficacies = make_gaussian_pmf(0, 1.5, 3)
        super (efficacies)
        # update based on an exam score
        update(score)
    end

    def likelihood(data, hypo)
        # Computes the likelihood of a test score, given efficacy.
        efficacy = hypo
        score = data
        raw = @exam.reverse(score)

        pmf = @exam.pmf_correct(efficacy)
        like = pmf.prob(raw)
        return like
    end

    def make_predictive_dist()
        # Returns the distribution of raw scores expected on a re-test.
        raw_pmf = @exam.make_raw_score_dist(self)
        return raw_pmf
    end
    
    def plot_posteriors(other)
        # Plots posterior distributions of efficacy.
		#    self, other: Sat objects.
        cdf1 = make_cdf_from_pmf(self)
        cdf2 = make_cdf_from_pmf(other)
		d1 = cdf1.items.transpose
		d2 = cdf2.items.transpose
		Numo.gnuplot do
#			set axis: "0.7, 1.0, 0.0, 1.0"
       		set title: 'sat_posteriors_eff'
			set xlabel: 'efficacy'
			set ylabel: "cdf"
		    plot [d1[0], d1[1], w:'lines', t:'self'], [d2[0], d2[1], w:'lines', t:'other']
		    sleep (2)
    	end
	end
end

def plot_joint_dist(pmf1, pmf2, thresh = 0.8)
    # Plot the joint distribution of p_correct.
	#    pmf1, pmf2: posterior distributions
    #    thresh: lower bound of the range to be plotted
    def clean(pmf, thresh)
        # Removes values below thresh.
        vals = pmf.values.select{|val| val < thresh }
        vals.each{|val| pmf.remove(val) }
	end
    clean(pmf1, thresh)
    clean(pmf2, thresh)
    pmf = make_joint(pmf1, pmf2)
	make_plot_3d(pmf.render, x_label: 'p_correct Alice', y_label: 'p_correct Bob',
					 title: 'sat_joint', contour: "both")
end

def compare_posterior_predictive(a_sat, b_sat)
    # Compares the predictive distributions of raw scores.
	#   a_sat: posterior distribution
    a_pred = a_sat.make_predictive_dist()
    b_pred = b_sat.make_predictive_dist()

	make_plot_2d(a_pred.render, x_label: "sat_score", y_label: "PMF", title: "prior")
	make_plot_2d(b_pred.render, x_label: "sat_score", y_label: "PMF", title: "posterior")

    print "Posterior predictive\n"
    print 'A = ', pmf_prob_greater(a_pred, b_pred), "\n"
    print 'B = ', pmf_prob_less(a_pred, b_pred), "\n"
    print 'C = ', pmf_prob_equal(a_pred, b_pred), "\n"
end

def plot_prior_dist(pmf)
    cdf = make_cdf_from_pmf(pmf)
	make_plot_2d(cdf.render, x_label: 'p_correct', y_label: "CDF", title: "prior")
end

class TopLevel < PmfSuite
    # Evaluates the top-level hypotheses about Alice and Bob.
	# Uses the bottom-level posterior distribution about p_correct (or efficacy).
	def initialize(value = 'None')
		super(value)
	end

    def update(data)
        a_sat, b_sat = data

        a_like = pmf_prob_greater(a_sat, b_sat)
        b_like = pmf_prob_less(a_sat, b_sat)
        c_like = pmf_prob_equal(a_sat, b_sat)

        a_like += c_like / 2
        b_like += c_like / 2

        mult('A', a_like)
        mult('B', b_like)

        self.normalize()
	end
end

def prob_correct(efficacy, difficulty, a = 1)
    # Returns the probability that a person gets a question right.
	#   efficacy: personal ability to answer questions
    #   difficulty: how hard the question is
	#   Returns: float prob"
    return 1 / (1 + Math.exp(-a * (efficacy - difficulty)))
end

def binary_pmf(p)
    # Makes a Pmf with values 1 and 0.
    #   p: probability given to 1
    #   Returns: Pmf object
    pmf = PmfSuite.new()
    pmf.set(1, p)
    pmf.set(0, 1 - p)
    return pmf
end

def pmf_correct(efficacy, difficulties)
    # Computes the distribution of correct responses.
	#   efficacy: personal ability to answer questions
    #   difficulties: list of difficulties, one for each question
	#   Returns: new Pmf object
    pmf0 = PmfSuite.new([0])
    ps = difficulties.collect{|difficulty| prob_correct(efficacy, difficulty) }
    pmfs = ps.collect{|p| binary_pmf(p) }
    dist = pmfs.inject(pmf0){|a, b| a + b }		#sum(pmfs, pmf0)
    return dist
end
alias toplevel_pmf_correct  pmf_correct

def make_difficulties(center, width, n)
    # Makes a list of n difficulties with a given center and width.
	#    Returns: list of n floats between center-width and center+width
    low, high = center - width, center + width
    return low.step(high, width * 2.0 / (n - 1)).to_a
end

def prob_correct_table()
    # Makes a table of p_correct for a range of efficacy and difficulty.
    efficacies = [3, 1.5, 0, -1.5, -3]
    difficulties = [-1.85, -0.05, 1.75]

    efficacies.each{|eff| 
        print '%5.2f || ' % eff
        difficulties.each{|diff|
            pr = prob_correct(eff, diff)
            print '%0.2f & ' % pr
        }
        print "\n"
    }
end
############################################
##  メインルーチン
############################################
prob_correct_table()
exam = Exam.new()
plot_prior_dist(exam.prior)
exam.calibrate_difficulty()
exam.compare_scores(780, 740, constructor = "Sat")
exam.compare_scores(780, 740, constructor = "Sat2")
