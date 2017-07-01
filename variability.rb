#! ruby -Ks
require './thinkbayes'
require "numo/gnuplot"
require 'date'

NUM_SIGMAS = 1

class Height < Joint
	# Hypotheses about parameters(mu,sigma) of the distribution of height.
    def initialize(mus, sigmas, name: nil)
		#    mus: sequence of possible mus
        # sigmas: sequence of possible sigmas
        pairs = []
		mus.each{|num| sigmas.each{|sigma| pairs << [num, sigma]}}.to_a
        super(pairs, name: name)
  	end

    def likelihood(data, hypo)
        mu, sigma = hypo
        like = eval_gaussian_pdf(data, mu, sigma)
        return like
	end

    def log_likelihood(data, hypo)
        mu, sigma = hypo
        return eval_gaussian_log_pdf(data, mu, sigma)
	end

    def log_update_set_fast(xs)
		# Updates the suite using a faster implementation.
		# Computes the sum of the log likelihoods directly.
        n = xs.size
        self.values().each{|hypo|
            mu, sigma = hypo
            total = summation(xs, mu)
            loglike = -n * Math.log(sigma) - total / 2 / sigma**2
            incr(hypo, loglike)
		}
	end

    def log_update_set_mean_var(data)
		#   data: sequence of values
        log_update_set_ABC(data.size, data.mean, data.sd)
	end

    def log_update_set_median_IPR(data)
		#   data: sequence of values
        median, s = median_s(data, num_sigmas = NUM_SIGMAS)
        print "median, s = ", median, ", ", s, "\n"
        log_update_set_ABC(data.size, median, s)
	end

    def log_update_set_ABC(n, m, s)
		#   n: sample size
        #   m: estimated central tendency
        #   s: estimated spread
        values.sort.each{|hypo|
            mu, sigma = hypo
            # compute log likelihood of m, given hypo
            stderr_m = sigma / Math.sqrt(n)
            loglike = eval_gaussian_log_pdf(m, mu, stderr_m)
            #compute log likelihood of s, given hypo
            stderr_s = sigma / Math.sqrt(2 * (n-1))
            loglike += eval_gaussian_log_pdf(s, sigma, stderr_s)

            incr(hypo, loglike)
		}
	end
end

def eval_gaussian_log_pdf(x, mu, sigma)
	#  x: float values
    #  mu, sigma: paramemters of Gaussian
    return Math.log(eval_gaussian_pdf(x, mu, sigma))
end

def find_prior_ranges(xs, num_points, median_flag = false)
	num_stderrs = 3.0
	# mu と sigma の尤度走査範囲を決める。
	#  xs: 与えれたsampleデータの分布を元に計算する
    #  num_points: データ数
    #  num_stderrs: 操作範囲(標準偏差の倍率)
    n = xs.size
    if median_flag
        m, s = median_s(xs, num_sigmas = NUM_SIGMAS)
    else
        m = xs.mean
        s = xs.sd
	end
    print 'classical estimators = ', m, "  ", s, "\n"

    # compute ranges for m(平均) and s(標準偏差)
    stderr_m = s / Math.sqrt(n)
    spread = stderr_m * num_stderrs
	step = spread * 2.0 / (num_points - 1)
    mus = (m - spread).step(m + spread, step).to_a

    stderr_s = s / Math.sqrt(2 * (n-1))
    spread = stderr_s * num_stderrs
	step = spread * 2.0 / (num_points - 1)
    sigmas = (s - spread).step(s + spread, step).to_a

    return mus, sigmas
end

def summation(xs, mu)
    # Computes the sum of (x-mu)**2 for x in t.
    ds = xs.collect{|x| (x - mu) ** 2}
    total = ds.sum
    return total
end

def coef_variation(suite)
    # 変動係数(標準偏差/平均)の確率マップ作成
	#     suite: Pmf that maps (x, y) to z
    pmf = PmfSuite.new()
    suite.items.each{|x, p| m, s = x
        pmf.incr(s / m, p) }
    return pmf
end

def plot_posterior(suite, pcolor = false, contour = true)
	dd = suite.collect{|v| v.flatten}.transpose
	Numo.gnuplot do
		set title:"Posterior joint distribution"
		set xlabel:"Mean height (cm)"
  		set ylabel:"Stddev (cm)"
		splot dd[0], dd[1], dd[2]#, w:'pm3d'
	sleep (2)
	end
end

def plot_coef_variation(suites)
    # 変動係数についての事後確率分布プロット
	t = []
    pmfs = {}
    suites.each{|label, suite|
        pmf = coef_variation(suite)
        print '#{label} CV posterior mean = ', pmf.mean(), "\n"
        cdf = make_cdf_from_pmf(pmf)
        pmfs[label] = pmf
      	t << [cdf.xs, cdf.ps, {w:'lines', t:label}]
	}
	make_plot_2d(t, x_label: 'Coefficient of variation',
					y_label: 'Probability', title: 'variability_cv', grid: "on", )
    print 'female bigger = ', pmf_prob_greater(pmfs['female'], pmfs['male']), "\n"
    print 'male bigger = ', pmf_prob_greater(pmfs['male'],pmfs['female']), "\n"
end

def plot_marginals(suite)
    # Plots marginal distributions from a joint distribution.
	#  suite: joint distribution of mu and sigma.
    cdf_m = make_cdf_from_pmf(suite.marginal(0))
    cdf_s = make_cdf_from_pmf(suite.marginal(1))

    Numo.gnuplot do
    	set title: '平均の分布'
		set xlabel: 'xy位置'
		set ylabel: 'CDF'
		plot cdf_m.xs, cdf_m.ps, {w:'lines', t:cdf_m.name}
	    sleep (5)
    end
    Numo.gnuplot do
    	set title: '標準偏差の分布'
		plot cdf_s.xs, cdf_s.ps, {w:'lines', t:cdf_s.name}
	    sleep (5)
    end
end

##################   以下ベイス更新ルーチン郡    ###################
#    suite: Suite that maps from (mu, sigma) to prob
#       xs: sequence
def update_suite1(suite, xs)
    suite.update_set(xs)
end

def update_suite2(suite, xs)
    suite.log.log_update_set(xs)
    suite.exp.normalize
end

def update_suite3(suite, xs)
    suite.log.log_update_set_fast(xs)
    suite.exp.normalize
end

def update_suite4(suite, xs)
    suite.log.log_update_set_mean_var(xs)
    suite.exp.normalize
end

def update_suite5(suite, xs)
    suite.log.log_update_set_median_IPR(xs)
    suite.exp.normalize
end

def median_IPR(xs, p)
    # Computes the median and interpercentile range.
	#   xs: sequence of values
    #    p: range (0-1), 0.5 yields the interquartile range
	#  returns: tuple of float (median, IPR)
    cdf = make_cdf_from_list(xs)
    median = cdf.percentile(50)
    alpha = (1-p) / 2
    ipr = cdf.value(1-alpha) - cdf.value(alpha)
    return median, ipr
end

def median_s(xs, num_sigmas)
    # Computes the median and an estimate of sigma.
	# Based on an interpercentile range (IPR).
    half_p = standard_gaussian_cdf(num_sigmas) - 0.5
    median, ipr = median_IPR(xs, half_p * 2)
    s = ipr / 2 / num_sigmas
    return median, s
end

def summarize(xs)
    xs.sort!
    print 'smallest = ', xs[0, 9], "\n"
    print 'largest = ', xs[-9, 10], "\n"
    cdf = make_cdf_from_list(xs)
    print "25%=", cdf.percentile(25), ", 50%=", cdf.percentile(50), ", 75%=", cdf.percentile(75), "\n"
end

def run_estimate(update_func, num_points = 31, median_flag = false)
    # Runs the whole analysis.
	#     update_func: which of the update functions to use
    #      num_points: number of points in the Suite (in each dimension)
    # DL出来ないので、データを製作。　男 平均:178.0、標準偏差:7.7、　女　平均:163.2、標準偏差:7.3
    d = [[1,[]],[2,[]]]
    rd = Random.new
    50000.times{ d[0][1] << rd.normal(178.0, 7.7) }
    50000.times{ d[1][1] << rd.normal(163.2, 7.3) }
    labels = {1 => 'male', 2 => 'female'}

    suites = {}
    d.each{|key, xs|
        name = labels[key]
        print name, "  ", xs.size, "\n"
        summarize(xs)

        xs = jitter(xs, 1.3)

        mus, sigmas = find_prior_ranges(xs, num_points, median_flag)
        suite = Height.new(mus, sigmas)
        suites[name] = suite
        send(update_func, suite, xs)
        print 'MLE = ', suite.maximum_likelihood(), "\n"

        plot_posterior(suite)

        pmf_m = suite.marginal(0)
        pmf_s = suite.marginal(1)
        print 'marginal mu_ave = ', pmf_m.mean(), 'mu_var = ', pmf_m.var(), "\n"
        print 'marginal sigma_ave = ', pmf_s.mean(), 'sigma_var = ', pmf_s.var(), "\n"
        plot_marginals(suite)
	}
    plot_coef_variation(suites)
end

##############################################################
###  mainルーチン
##############################################################
func = :update_suite5
#func = :update_suite4
#func = :update_suite3
#func = :update_suite2
#func = :update_suite1

ini = Time.now
median_flag = (func == :UpdateSuite5)
run_estimate(func, 31, median_flag)
print "#{func}　所要時間 =", Time.now - ini, "\n"
