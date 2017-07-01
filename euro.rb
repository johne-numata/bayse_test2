
puts "Hello World!"

require "./Thinkbayes.rb"
require "./Thinkbayes_2.rb"

class Euro < Suite

	def likelihood(data, hypo)
		x = hypo
        if data == "H" || data == "h"
        	return x / 100.0
        else
        	return 1.0 - x / 100.0
        end
	end
end

hypos = (0..100).to_a
suite = Euro.new(hypos)

#suite.print_dict

dataset = ["h"] * 140 + ["t"] * 110
dataset.each{|x| suite.update(x)}

suite.print_dict

p suite.maximumLikelihood
p suite.mean
p percentile(suite, 50)
p credibleInterval(suite, 90)
p suite.prob(50)
