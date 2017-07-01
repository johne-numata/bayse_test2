#! ruby -Ks
require './thinkbayes'
require "numo/gnuplot"

class Emitter < PmfSuite
	# Represents hypotheses about r.
    def initialize(rs, f = 0.1)
        # rs: sequence of hypothetical emission rates
        # f: fraction of particles registered
        detectors = rs.collect{|r| Detector.new(r, f) }
        super(detectors)
	end

    def update(data)
        # Updates the Suite based on data.
		#    data: number of particles counted
        super(data)  # 親クラスのupdate呼び出し
        self.values.each {|detector| detector.update(data) }
    end

    def likelihood(data, hypo)
        # Likelihood of the data given the hypothesis.
		#   data: number of particles counted
        #   hypo: emission rate, r
		#   Returns: probability density of the data under the hypothesis
        detector = hypo
        like = detector.suite_likelihood(data)
        return like
	end

    def dist_of_r()
        # Returns the PMF of r.
        items = self.items.collect{|detector, prob| [detector.r, prob] }
        return make_suite_from_items(items)
	end

    def dist_of_n()
        # Returns the PMF of n.        
        return make_mixture(self)
    end
end

class Emitter2 < PmfSuite
    # Represents hypotheses about r.
	def initialize(rs, f = 0.1)
        # rs: sequence of hypothetical emission rates
        # f: fraction of particles registered
        detectors = rs.collect{|r| Detector.new(r, f) }
        super(detectors)
	end

    def likelihood(data, hypo)
        # Likelihood of the data given the hypothesis.
		#   data: number of counted per unit time
        #   hypo: emission rate, r
		#   Returns: probability density of the data under the hypothesis
		t = hypo.update(data)
        return t
#        return hypo.update(data)
    end

    def dist_of_r()
        # Returns the PMF of r.
        items = self.items.collect{|detector, prob| [detector.r, prob] }
        return make_suite_from_items(items)
    end

    def dist_of_n()
        # Returns the PMF of n.        
        return make_mixture(self)
	end
end

class Detector < PmfSuite
    # Represents hypotheses about n.
    attr_reader  :r
	def initialize(r, f, high = 500, step = 5)
		# r: known emission rate, r
        # f: fraction of particles registered
        # high: maximum number of particles, n
        # step: step size between hypothetical values of n
        pmf = make_poisson_pmf(r, high, step)
        super(pmf)
        @r = r
        @f = f
	end

    def likelihood(data, hypo)
        # Likelihood of the data given the hypothesis.
		#   data: number of particles counted
        #   hypo: number of particles hitting the counter, n
        k = data
        n = hypo
        p = @f
        return eval_binomial_pmf(k, n, p)
	end

    def suite_likelihood(data)
        # Adds up the total probability of the data under the suite.
		#   data: number of particles counted
        total = 0
        self.items.each{|hypo, prob|
            like = likelihood(data, hypo)
            total += prob * like
        }
        return total
    end
end
##########################################################
###   以下mainルーチン
##########################################################
k = 15
f = 0.1

# plot Detector suites for a range of hypothetical r
t = [100,200,400].collect{|r| 
	suite = Detector.new(r, f, 500, 1)
    suite.update(k)
    print suite.maximum_likelihood(), "\n"
    suite.render
    }
make_plot_2d(t, x_label: 'Number of particles (n)', y_label: 'PMF', title: 'jaynes1', grid: "on", )

# plot the posterior distributions of r and n
hypos = 1.step(501, 5).to_a
suite = Emitter2.new(hypos, f)
suite.update(k)

post_r = suite.dist_of_r(); post_r.name = "post_r"
post_n = suite.dist_of_n(); post_n.name = "post_n"

make_plot_2d([post_r.render, post_n.render], x_label: 'Emission rate', y_label: "pmf", title: 'jaynes2', grid: "on", )
