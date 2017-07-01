
puts "Hello World!"

require "./Thinkbayes.rb"

class Monty < Pmf
	def initialize(hypos = "None")
		super
	    return  if hypos == "None"
        hypos.each{|x| set(x, 1)}
        normalize
    end

	def likelihood(data, hypo)
        if hypo == data
            return 0
        elsif hypo == 'A'
            return 0.5
        else
            return 1
        end
	end

	def update(data)
        @d.each_key{|hypo|
            like = self.likelihood(data, hypo)
            self.mult(hypo, like)
        }
        self.normalize()
	end

end

hypos = ['A', 'B', 'C']
pmf = Monty.new(hypos)
pmf.update('B')
p pmf
pmf.print_dict


#pmf = Pmf.new()
#words.each{|x| pmf.incr(x, 1) }
#pmf.set('Bowl 1', 0.5)
#pmf.set('Bowl 2', 0.5)
#p pmf
#pmf.mult('Bowl 1', 0.75)
#pmf.mult('Bowl 2', 0.5)
#pmf.normalize
#p pmf
#p pmf.prob("this")
#p pmf.maxlike