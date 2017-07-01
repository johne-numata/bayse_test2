
puts "Hello World!"

require "./Thinkbayes.rb"

class M_and_M < Suite
	def initialize (value = "None")
		super
		mix94 = {"brown" => 30, "yellow" => 20, "red" => 20,
                 "green" => 10, "orange" => 10, "tan" => 10}
	    mix96 = {"blue" => 24, "green" => 20, "orange" => 16,
                 "yellow" => 14, "red" => 13, "brown" => 13}

    	hypoA = {"bag1" => mix94, "bag2" => mix96}
    	hypoB = {"bag1" => mix96, "bag2" => mix94}

    	@hypotheses = {"A" => hypoA, "B" => hypoB}
	end
	
	def likelihood(data, hypo)
        	### Computes the likelihood of the data under the hypothesis.
			### hypo: string hypothesis (A or B)
        	### data: tuple of string bag, string color
        bag, color = data
        mix = @hypotheses[hypo][bag]
        like = mix[color]
        return like
	end

end

hypos = ['A', 'B']
suite = M_and_M.new(hypos)

suite.update(['bag1', 'yellow'])
suite.print_dict
suite.update(['bag2', 'green'])
suite.print_dict

