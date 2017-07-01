#! ruby -Ks
require './thinkbayes'
require "numo/gnuplot"

def strafing_speed(alpha, beta, x)
    # 掃射速度の計算、角度が浅いと着弾点の移動速度が早い
	# alpha: x location of shooter
    # beta: y location of shooter
    # x: location of impact
	# Returns: derivative of x with respect to theta
    theta = Math.atan2(x.to_f - alpha.to_f, beta.to_f)
    speed = beta.to_f / Math.cos(theta.to_f)**2
    return speed
end

def make_location_pmf(alpha, beta, locations)
	# alpha/betaの地点にいた場合の、着弾点別の確率をPmfで表す。
	# 確率は掃射速度に反比例するとして計算する。
	# alpha: x position
    # beta: y position
    # locations: x locations where the pmf is evaluated
    pmf = PmfSuite.new(name: beta.to_s)
    locations.each{|x|
        prob = 1.0 / strafing_speed(alpha, beta, x)
        pmf.set(x, prob)
    }
    pmf.normalize()
    return pmf
end

class Paintball < Joint
    # 相手の位置に関する仮説を表します
    def initialize(alphas, betas, locations)
        # Stores locations for use in Likelihood.
		#   alphas: 存在することが可能なx位置の配列
        #   betas: 存在することが可能なy位置の配列
        #   locations:　着弾可能な位置の配列
		@pairs = []
        @locations = locations
        alphas.each{|a| betas.each{|b| @pairs << [a, b]}}
        super (@pairs)
	end

    def likelihood(data, hypo)
		# Computes the likelihood of the data under the hypothesis.
		#   hypo: pair of alpha, beta
        #   data: location of a hit
		#   Returns: float likelihood
        alpha, beta = hypo
        x = data
        pmf = make_location_pmf(alpha, beta, @locations)
        like = pmf.prob(x)
        return like
	end
end

def make_pmf_plot(alpha = 10)
    # Plots Pmf of location for a range of betas.
    locations = (0..31).to_a
    betas = [10, 20, 40]
	pmfs = []
    betas.each{|beta| pmfs << make_location_pmf(alpha, beta, locations) }

	t = []
	pmfs.each{|pmf|
		d = pmf.items.transpose
      	t << [d[0], d[1], {w:'lines', t: pmf.name}]
	}
	make_plot_2d(t, x_label: "location", y_label: "pmf", title: "beta依存",	grid: "on", )
end

def make_posterior_plot(suite)
    # Plots the posterior marginal distributions for alpha and beta.
    marginal_alpha = suite.marginal(0)
    marginal_beta = suite.marginal(1)
    print 'alpha CI ', marginal_alpha.credible_interval(50), "\n"
    print 'beta CI ', marginal_beta.credible_interval(50), "\n"

	t = [marginal_alpha.render, marginal_beta.render]
	make_plot_2d(t, x_label: "location", y_label: "pmf", title: "事後確率", grid: "on")

	t = [marginal_alpha.make_cdf.render, marginal_beta.make_cdf.render]
	make_plot_2d(t, x_label: "location", y_label: "cdf", title: "事後確率", grid: "on")
end

def make_conditional_plot(suite)
    # Plots marginal CDFs for alpha conditioned on beta.
	# suite: posterior joint distribution of location
    betas = [10, 20, 40]
	t = []
    betas.each{|beta| t << suite.conditional(0, 1, beta).render }
	make_plot_2d(t, x_label: "location", y_label: "pmf", title: "beta条件別alpha積算した事後確率", grid: "on", )
end

def make_contour_plot(suite)
    # Plots the posterior joint distribution as a contour plot.
	dd = suite.items.collect{|v| v.flatten}.transpose
   	t = [dd[0], dd[1], dd[2], {w:'pm3d at s', t: suite.name}]
	make_plot_3d(t, x_label: "alpha", y_label: "beta", title: "paintball_14",
					xrange: "[0:30]", yrange: "[0:20]", contour: "both")
end

def make_credible_plot(suite)
    # Makes a plot showing several two-dimensional credible intervals.
    d = suite.values.collect{|pair| [pair, 0]}.to_h
    percentages = [75, 50, 25]
    percentages.each{|p|
        interval = suite.max_like_interval(p)
        interval.each{|pair|
        	pair
            d[pair] += 1
        }
	}

	dd = d.collect{|v| v.flatten}.transpose
   	t = [dd[0], dd[1], dd[2], {w:'pm3d'}]
	make_plot_3d(t, x_label: "alpha", y_label: "beta", title: "確率Map", grid: "on", map: "on",
					xrange: "[0:31]", yrange: "[0:51]")
end

###############################################################
###    以下 mainスクリプト
###############################################################
alphas = (0..31).to_a
betas = (1..51).to_a
locations = (0..31).to_a

suite = Paintball.new(alphas, betas, locations)
suite.update_set([15, 16, 18, 21])

make_credible_plot(suite)
make_contour_plot(suite)
make_posterior_plot(suite)
make_conditional_plot(suite)
make_pmf_plot()
