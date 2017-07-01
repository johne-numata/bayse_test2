#! ruby -Ks

require 'csv'
require "./kde.rb"
require "./thinkbayes"
require "numo/gnuplot"

class Price < PmfSuite
    # Represents hypotheses about the price of a showcase.
    def initialize(pmf, player)
		# pmf: prior distribution of price
        # player: Player object
        super(pmf)
        @player = player
    end

    def likelihood(data, hypo)
        # Computes the likelihood of the data under the hypothesis.
		#    hypo: actual price
        #    data: the contestant's guess
	    price = hypo
        guess = data

        error = price - guess
        like = @player.error_density(error)

        return like
	end
end

class GainCalculator
    # Encapsulates computation of expected gain.
    def initialize(player, opponent)
        @player = player
        @opponent = opponent
	end

    def expected_gains(bid_xs = nil)
		# Computes expected gains for a range of bids.
    	#    returns: tuple (sequence of bids, sequence of gains)
		gains = []
		bid_xs = @player.price_xs  if bid_xs.nil?
		bid_xs.each {|bid| gains << expected_gain(bid) }
        return bid_xs, gains
	end

    def expected_gain(bid)
        suite = @player.posterior
        total = 0
        for price, prob in suite.items().sort
            gain_a = gain(bid, price)
            total += prob * gain_a
        end
        return total
    end

    def gain(bid, price)
    	# Computes the return of a bid, given the actual price.
		#    bid: number
        # 	 price: actual price

        # if you overbid, you get nothing
        return 0 if bid > price

        # otherwise compute the probability of winning
        diff = price - bid
        prob = prob_win(diff)

        # if you are within 250 dollars, you win both showcases
        if diff <= 250
            return 2 * price * prob
        else
            return price * prob
        end
	end

    def prob_win(diff)
        # Computes the probability of winning for a given diff.
		#    diff: how much your bid was off by
        prob = (@opponent.prob_overbid() + 
                @opponent.prob_worse_than(diff))
        return prob
	end
end

class Player
    # Represents a player on The Price is Right.
	attr_reader :prior, :posterior, :price_xs, :cdf_diff

    def initialize(prices, bids, diffs)
        @pdf_price = EstimatedPdf.new(prices)
        @cdf_diff = make_cdf_from_list(diffs)
	    @price_xs = 0.step(75000, 750).to_a   # 101個のデータを生成
	    
        mu = 0
        sigma = diffs.sd
        @pdf_error = GaussianPdf.new(mu, sigma)
	end

    def error_density(error)
        # Density of the given error in the distribution of error.
        # error: how much the bid is under the actual price
        return @pdf_error.density(error)
	end

    def pmf_price()
        # Returns a new Pmf of prices.
		# A discrete version of the estimated Pdf.
        return @pdf_price.make_pmf(@price_xs)
	end

    def prob_overbid()
        # Returns the probability this player overbids.
        return @cdf_diff.prob(-1)
	end

    def prob_worse_than(diff)
        # Probability this player's diff is greater than the given diff.
		# diff: how much the oppenent is off by (always positive)
        return 1 - @cdf_diff.prob(diff)
	end

    def make_beliefs(guess)
        # Makes a posterior distribution based on estimated price.
		# Sets attributes prior and posterior.
		# guess: what the player thinks the showcase is worth        
        pmf = pmf_price()
        @prior = Price.new(pmf, self)
        @posterior = Marshal.load(Marshal.dump(@prior))  #オブジェクトの深いコピー
        @posterior.update(guess)
	end

    def optimal_bid(guess, opponent)
        # Computes the bid that maximizes expected return.
        #   guess: what the player thinks the showcase is worth 
        #   opponent: Player
		#   Returns: (optimal bid, expected gain)
        make_beliefs(guess)
        calc = GainCalculator.new(self, opponent)
        bids, gains = calc.expected_gains()
        gain, bid = gains.zip(bids).max{|a, b| a[0] <=> b[0]}
        return bid, gain
    end

    def plot_beliefs(root)
		# Plots prior and posterior beliefs.
		# 	root: string filename root for saved figure
	    pre = @prior.items().transpose()
    	t = [pre[0], pre[1], w:'lines', t: 'prior']
	    post = @posterior.items().transpose()
    	t << [post[0], post[1], w:'lines', t: 'posterior']
	    post_name = 'posterior'
		make_plot_2d(t, x_label: "price($)", y_label: "PDF", title: "事前/事後予測", grid: "on", )
	end			
end

def make_plots(player_1, player_2)
    # price1 shows the priors for the two players
    # price2 shows the distribution of diff for the two players

    # plot the prior distribution of price for both players
    pmf1 = player_1.pmf_price().items().transpose()
    t = [pmf1[0], pmf1[1], w:'lines', t: 'showcase 1']
    pmf2 = player_2.pmf_price().items().transpose()
    t << [pmf2[0], pmf2[1], w:'lines', t: 'showcase 2']
	make_plot_2d(t, x_label: "price($)", y_label: "PDF", title: "ショーケースの価格分布", grid: "on", )

    # plot the historical distribution of underness for both players
    print 'Player diff median = ', player_1.cdf_diff.percentile(50),"\n"
    print 'Player diff median = ', player_2.cdf_diff.percentile(50),"\n"
    print 'Player 1 overbids prob = ', player_1.prob_overbid(),"\n"
    print 'Player 2 overbids prob = ', player_2.prob_overbid(),"\n"
    t = [player_1.cdf_diff.xs, player_1.cdf_diff.ps, w:'lines', t: 'player 1']
    t << [player_2.cdf_diff.xs, player_2.cdf_diff.ps, w:'lines', t: 'player 2']
	make_plot_2d(t, x_label: "diff($)", y_label: "CDF", title: "推測値と実価格の差異の累積分布", grid: "on", )
	
end

def make_players()
    # Reads data and makes player objects.
  array = CSV.read("./showcases.2011.csv")
  array2 = CSV.read("./showcases.2012.csv")
  array.select!{|a| a.compact != []}
  array2.select!{|a| a.compact != []}
  array.each_index{|idx| array[idx] = array[idx] + array2[idx].drop(1)}

  price_1 = array[2][1..array[2].size-1]
  price_2 = array[3][1..array[2].size-1]
  bid_1 = array[4][1..array[2].size-1]
  bid_2 = array[5][1..array[2].size-1]
  diff_1 = array[6][1..array[2].size-1]
  diff_2 = array[7][1..array[2].size-1]

  player_1 = Player.new(price_1, bid_1, diff_1)
  player_2 = Player.new(price_2, bid_2, diff_2)

  return player_1, player_2
end

def plot_expected_gains(guess_1=20000, guess_2=40000)
    # Plots expected gains as a function of bid.
    # guess1: player1's estimate of the price of showcase 1
    # guess2: player2's estimate of the price of showcase 2
    
    player_1, player_2 = make_players()
    make_plots(player_1, player_2)

    player_1.make_beliefs(guess_1)
    player_2.make_beliefs(guess_2)

    print "Player 1 prior mle = ", player_1.prior.maximum_likelihood(), "\n"
    print "Player 2 prior mle = ", player_2.prior.maximum_likelihood(), "\n"
    print "Player 1 posterior mean = ", player_1.posterior.mean(), "\n"
    print "Player 2 posterior mean = ", player_2.posterior.mean(), "\n"
    print "Player 1 posterior mle = ", player_1.posterior.maximum_likelihood(), "\n"
    print "Player 2 posterior mle = ", player_2.posterior.maximum_likelihood(), "\n"

    player_1.plot_beliefs('price3')
    player_2.plot_beliefs('price4')

	calc_1 = GainCalculator.new(player_1, player_2)
	calc_2 = GainCalculator.new(player_2, player_1)

    bids, gains = calc_1.expected_gains()
    t = [bids, gains, w:'lines', t: 'player 1']
    print 'Player 1 optimal bid = ', gains.zip(bids).max{|a, b| a[0] <=> b[0]}, "\n"
    bids2, gains2 = calc_2.expected_gains()
    t << [bids2, gains2, w:'lines', t: 'player 2']
    print 'Player 2 optimal bid = ', gains2.zip(bids2).max{|a, b| a[0] <=> b[0]}, "\n"
	make_plot_2d(t, x_label: "bid ($)", y_label: "expected gain ($)",
							 title: "期待値の分布", grid: "on", )
end

def plot_optimal_bid()
    # Plots optimal bid vs estimated price.
    player_1, player_2 = make_players()
    guesses = 15000.step(120000, 2250).to_a
	temp = Marshal.dump(player_1)		 #オブジェクトの深いコピー用
	
    res = []
    for guess in guesses
    	player = Marshal.load(temp) 
		bid, gain = player.optimal_bid(guess, player_2)
        mean = player.posterior.mean()
        mle = player.posterior.maximum_likelihood()
        res << [guess, mean, mle, gain, bid]
    end

    labels = ["guess", "mean", "mle", "gain", "bid"]
	res = res.transpose
	t = []
	res.each_index{|i| t << [guesses, res[i], w:'lines', t: labels[i]] if i > 0}
	make_plot_2d(t, x_label: "guessed price ($)", title: "予測値別期待値", grid: "on", file_name: "test.gif")
end

def test_code(calc)
    # Check some intermediate results.
	#    calc: GainCalculator

    # test ProbWin
    for diff in [0, 100, 1000, 10000, 20000]
        print diff, calc.prob_win(diff)
    end

    # test Return
    price = 20000
    for bid in [17000, 18000, 19000, 19500, 19800, 20001]
        print bid, calc.gain(bid, price)
    end
end

plot_expected_gains()
plot_optimal_bid()
