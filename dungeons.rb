#! ruby -Ks
require "./Thinkbayes.rb"

class Die < PmfSuite
	def initialize(sides = 1.0)
		super()
		sides.times{|x| self.set(x + 1, 1)}
		self.normalize()
	end

end

def pmf_max(pmf1, pmf2)	# Computes the distribution of the max of values drawn from two Pmfs.
						# pmf1, pmf2: Pmf objects
						# returns: new Pmf
    res = PmfSuite.new()
    pmf1.each{|v1, p1|
    	pmf2.each{|v2, p2|
            res.incr([v1, v2].max, p1 * p2)
		}
	}
	res.normalize
    return res
end

# ˆÈ‰ºMIX‚ÌŽÀK

 pmf_dice = PmfSuite.new()
 pmf_dice.set(Die.new(4), 5)
 pmf_dice.set(Die.new(6), 4)
 pmf_dice.set(Die.new(8), 3)
 pmf_dice.set(Die.new(12), 2)
 pmf_dice.set(Die.new(20), 1)
 pmf_dice.normalize()

# mix = PmfSuite.new()
# pmf_dice.each{|die, weight|
# 	 die.each{|outcome, prob|
# 	 	 mix.incr(outcome, weight * prob)
# 	} 
# }
 mix = make_mixture(pmf_dice)
 mix.print_dict

print "\n\n"
# ˆÈ‰ºMAX‚ÌŽÀK

 d6 = Die.new(6)
 dice = [d6] * 3
 three = sample_sum(dice, 10000)
 three.print_dict
print "\n"
 
 d6 = Die.new(6)
 dice = [d6] * 3
 three = sample_max(dice, 10000)
 three.print_dict
print "\n"

 three_exact = d6 + d6 + d6
 three_exact.print_dict
print "\n\n"

 best_attr2 = pmf_max(three_exact, three_exact)
 best_attr2.print_dict
 best_attr4 = pmf_max(best_attr2, best_attr2)
 best_attr4.print_dict
 best_attr6 = pmf_max(best_attr4, best_attr2)
 best_attr6.print_dict
print "\n\n"

print "\n\n"
 best_attr_cdf = three_exact.max(6)
 best_attr_pmf = best_attr_cdf.make_pmf
 best_attr_pmf.print_dict

