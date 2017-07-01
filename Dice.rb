
puts "Hello World!"

require "./Thinkbayes.rb"

class Dice < Suite
#	def initialize (value = "None")
#		super
#	end
	
	def likelihood(data, hypo)
        if hypo < data
        	return 0
        else
        	return 1.0 / hypo
        end
	end

end

hypos = [4, 6, 8, 12, 20]
suite = Dice.new(hypos)
suite.print_dict

suite.update(6)
suite.print_dict
[6, 8, 7, 7, 5, 4].each{|x| suite.update(x)} 
suite.print_dict

#suite.update(6)
#suite.print_dict

#suite.update(8)
#suite.print_dict

