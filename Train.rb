
puts "Hello World!"

require "./Thinkbayes.rb"
require "./Thinkbayes_2.rb"

class Train < Suite
	def initialize(values = "None", alpha = 1.0)
		super(values)
        values.each{|v| self.set(v, v**(-alpha))}
	end
	def likelihood(data, hypo)
        if hypo < data
        	return 0
        else
        	return 1.0 / hypo
        end
	end
end

hypos = (1..1000).to_a
suite = Train.new(hypos)

#suite.print_dict

suite.update(60)
suite.update(30)
suite.update(90)
suite.print_dict

p suite.mean
interval = percentile(suite, 5), percentile(suite, 95)
p interval

cdf = suite.makeCdf()
interval = cdf.percentile(5), cdf.percentile(95)
p interval
