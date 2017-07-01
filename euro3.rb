#! ruby -Ks
require './thinkbayes'
require "numo/gnuplot"

class Euro < PmfSuite
    # Represents hypotheses about the probability of heads.
	def likelihood(data, hypo)
		#   hypo: integer value of x, the probability of heads (0-100)
        #   data: tuple of (number of heads, number of tails)
        x = hypo / 100.0
        heads, tails = data
        like = x ** heads * (1-x) ** tails
        return like
	end
end

def triangle_prior()
    # Makes a Suite with a triangular prior.
    suite = Euro.new()
    0.step(51){|x| suite.set(x, x) }
    51.step(101){|x| suite.set(x, 100 - x) } 
    suite.normalize()
    return suite
end

def suite_likelihood(suite, data)
    # Computes the weighted average of likelihoods for sub-hypotheses.
	#    suite: Suite that maps sub-hypotheses to probability
    #     data: some representation of the data
    #  returns: float likelihood
    total = 0
    suite.items.each{|hypo, prob|
        like = suite.likelihood(data, hypo)
        total += prob * like
    }
    return total
end

############################################################
###  ƒƒCƒ“ƒ‹[ƒ`ƒ“
############################################################
data = 140, 110
#data = 8, 12

suite = Euro.new()
like_f = suite.likelihood(data, 50)
print 'p(D|F) = ', like_f, "\n"

actual_percent = 100.0 * 140 / 250
likelihood = suite.likelihood(data, actual_percent)
print 'p(D|B_cheat) = ', likelihood, "\n"
print 'p(D|B_cheat) / p(D|F = )', likelihood / like_f, "\n"
like40 = suite.likelihood(data, 40)
like60 = suite.likelihood(data, 60)
likelihood = 0.5 * like40 + 0.5 * like60
print 'p(D|B_two) = ', likelihood, "\n"
print 'p(D|B_two) / p(D|F) = ', likelihood / like_f, "\n"

b_uniform = Euro.new(0.step(100).to_a)
b_uniform.remove(50)
b_uniform.normalize()
likelihood = suite_likelihood(b_uniform, data)
print 'p(D|B_uniform) = ', likelihood, "\n"
print 'p(D|B_uniform) / p(D|F) = ', likelihood / like_f, "\n"

b_tri = triangle_prior()
b_tri.remove(50)
b_tri.normalize()
p likelihood = b_tri.update(data)
print 'p(D|B_tri) = ', likelihood, "\n"
print 'p(D|B_tri) / p(D|F) = ', likelihood / like_f, "\n"

