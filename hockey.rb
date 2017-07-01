#! ruby -Ks
require 'csv'
require './thinkbayes'
require "numo/gnuplot"

USE_SUMMARY_DATA = true

class Hockey < PmfSuite
	# Represents hypotheses about the scoring rate for a team."""

    def initialize()
        if USE_SUMMARY_DATA
            # prior based on each team's average goals scored
            mu = 2.8
            sigma = 0.3
        else
            # prior based on each pair-wise match-up
            mu = 2.8
            sigma = 0.85
		end
        pmf = make_gaussian_pmf(mu, sigma, 4)
        super(pmf)
    end
    
    def likelihood(data, hypo)
		# Computes the likelihood of the data under the hypothesis.
		# Evaluates the Poisson PMF for lambda and k.
		# hypo: goal scoring rate in goals per game
        # data: goals scored in one period
        lam = hypo
        k = data
        like = eval_poisson_pmf(k, lam)
        return like
    end
end

def make_goal_pmf(suite, high=10)
    # Makes the distribution of goals scored, given distribution of lam.
	# suite: distribution of goal-scoring rate
    # high: upper bound
	# returns: Pmf of goals per game

    metapmf = PmfSuite.new()

    suite.items.each{|lam, prob|
        pmf = make_poisson_pmf(lam, high)
        metapmf.set(pmf, prob)
	}
    mix = make_mixture(metapmf)
    return mix
end

def make_goal_time_pmf(suite)
    # Makes the distribution of time til first goal.
	# suite: distribution of goal-scoring rate
	# returns: Pmf of goals per game

    metapmf = PmfSuite.new()

    suite.items.each{|lam, prob|
        pmf = make_exponential_pmf(lam, high = 2, n = 2001)
        metapmf.set(pmf, prob)
	}
    mix = make_mixture(metapmf)
    return mix
end

def read_hockey_data(filename = 'hockey_data.csv')
	# Read game scores from the data file.
	# filename: string
    game_list = CSV.read(filename)

    # map from gameID to list of two games
    games = {}
    for game in game_list
        next if game.shift != "2011"
        key = game.shift
        games.fetch(key){|key| games[key] = []}
        games[key] << game
	end

    # map from (team1, team2) to (score1, score2)
    pairs = {}
    games.each{|key, pair|
        t1, t2 = pair
        key = t1[0], t2[0]
        entry = t1.last, t2.last
        pairs.fetch(key){|key| pairs[key] = []}
        pairs[key] << entry
	}
#p pairs	
    process_scores_teamwise(pairs)
    process_scores_pairwise(pairs)
end

def process_scores_pairwise(pairs)
    # Average number of goals for each team against each opponent.
	# pairs: map from (team1, team2) to (score1, score2)

    # map from (team1, team2) to list of goals scored
    goals_scored = {}
    pairs.each{|key, entries|
        t1, t2 = key
        key2 = t2, t1
        entries.each{|entry|
            g1, g2 = entry
            goals_scored.fetch(key){|key| goals_scored[key] = []}
	        goals_scored[key] << g1.to_i
            goals_scored.fetch(key2){|key| goals_scored[key] = []}
	        goals_scored[key2] << g2.to_i
        }
    }
#p goals_scored
    # make a list of average goals scored
    lams = []
    goals_scored.each{|key, goals|
		next if goals.size < 3
        lam = goals.mean
        lams << lam
	}
    # make the distribution of average goals scored
    cdf = make_cdf_from_list(lams).items.transpose
	Numo.gnuplot do
		set xlabel:"average goals"
  		set ylabel:"CDF"
		plot cdf[0], cdf[1], w:'lines'
	sleep (2)
	end

    mu, var = lams.mean, lams.var
    print 'mu, sig = ', mu, "  ", Math.sqrt(var), "\n"
    print 'BOS v VAN = ', pairs[['BOS', 'VAN']], "\n"
end

def process_scores_teamwise(pairs)
	# Average number of goals for each team.
	# pairs: map from (team1, team2) to (score1, score2)

    # map from team to list of goals scored
    goals_scored = {}
    pairs.each{|key, entries|
        t1, t2 = key
        entries.each{|entry|
            g1, g2 = entry
            goals_scored.fetch(t1){|key| goals_scored[key] = []}
	        goals_scored[t1] << g1.to_i
            goals_scored.fetch(t2){|key| goals_scored[key] = []}
	        goals_scored[t2] << g2.to_i
	    }
	}
#	p goals_scored
    # make a list of average goals scored
    lams = []
    goals_scored.each{|key, goals|
        lam = goals.mean
        lams << lam
	}
    # make the distribution of average goals scored
    cdf = make_cdf_from_list(lams).items.transpose
	Numo.gnuplot do
		set xlabel:"average goals"
  		set ylabel:"CDF"
		plot cdf[0], cdf[1], w:'lines'
	sleep (2)
	end

    mu, var = lams.mean, lams.var
    print 'mu, sig = ', mu, " ", Math.sqrt(var), "\n"
end

#########################################################################
###   以下メインルーチン
#########################################################################
read_hockey_data()
suite1 = Hockey.new()		# 'bruins'
suite2 = Hockey.new()		# 'canucks'

bruins = suite1.items().transpose()
canucks = suite2.items().transpose()
make_plot_2d([[bruins[0], bruins[1], w:'lines', t: 'bruins'], [canucks[0], canucks[1],
				 w:'lines', t: 'canucks']] , x_label: "goals per game",
				  y_label: "Probability", title: "事前分布", grid: "on", )

#過去4試合の結果でベイズ更新
suite1.update_set([0, 2, 8, 4])
suite2.update_set([1, 3, 1, 0])

bruins = suite1.items().transpose()
canucks = suite2.items().transpose()
make_plot_2d([[bruins[0], bruins[1], w:'lines', t: 'bruins'], [canucks[0], canucks[1],
				 w:'lines', t: 'canucks']] , x_label: "goals per game",
				  y_label: "Probability", title: "事後分布", grid: "on", )

#一試合のゴール数の計算
goal_dist1 = make_goal_pmf(suite1)
goal_dist2 = make_goal_pmf(suite2)

bruins = goal_dist1.items().transpose()
canucks = goal_dist2.items().transpose()
make_plot_2d([[bruins[0], bruins[1], w:'lines', t: 'bruins'], [canucks[0], canucks[1],
				 w:'lines', t: 'canucks']] , x_label: "goals",
				  y_label: "Probability", title: "次回のゴール獲得数の確率", grid: "on", )

#ゴール間隔の計算(何試合に1ゴールか)
time_dist1 = make_goal_time_pmf(suite1)    
time_dist2 = make_goal_time_pmf(suite2)
print 'MLE bruins = ', suite1.maximum_likelihood(), "\n"
print 'MLE canucks = ', suite2.maximum_likelihood(), "\n"
   
bruins = time_dist1.items().transpose()
canucks = time_dist2.items().transpose()
make_plot_2d([[bruins[0], bruins[1], w:'lines', t: 'bruins'], [canucks[0], canucks[1],
				 w:'lines', t: 'canucks']] , x_label: "Games until goal",
				  y_label: "Probability", title: "ゴールするまでの試合数の確率", grid: "on", )

#次の試合に勝つ確率/負ける確率の計算
diff = goal_dist1 - goal_dist2
p_win = diff.prob_greater(0)
p_loss = diff.prob_less(0)
p_tie = diff.prob(0)
print "p-win, p-loss, p-tie = ", p_win, "  ", p_loss, "  ",  p_tie, "\n"

#早くゴール出来る確率の計算、サドンデス検討用
p_overtime = pmf_prob_less(time_dist1, time_dist2)
p_adjust = pmf_prob_equal(time_dist1, time_dist2)
p_overtime += p_adjust / 2
print 'p_overtime = ', p_overtime, "\n" 
print "引き分けの内、サドンデスでかつ確率", p_overtime * p_tie, "\n"
p_win += p_overtime * p_tie
print '最終的な次の試合で勝つ確率　p_win = ', p_win, "\n"

# win the next two　2試合連続で勝つ確率
p_series = p_win**2

# split the next two, win the third　一試合負けて最終的に優勝の確率
p_series += 2 * p_win * (1-p_win) * p_win
print '現時点での最終的に優勝する確率　p_series = ', p_series, "\n"
