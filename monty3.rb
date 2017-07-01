
puts "Hello World!"

require "./Thinkbayes.rb"

class Monty < Suite
	def likelihood(data, hypo)
        if hypo == data
            return 0
        elsif hypo == 'A'
            return 0.5
        else
            return 1
        end
	end
end

hypos = ['A', 'B', 'C']
pmf = Monty.new(hypos)
p pmf
pmf.update('B')

pmf.print_dict

