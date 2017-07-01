
puts "Hello World!"

require "./Thinkbayes.rb"

class Cookie < Pmf
	def initialize(hypos = "None")
		super
	    return  if hypos == "None"
        hypos.each{|x| set(x, 1)}
        normalize
    end

    @@mixes = {
   	    'Bowl 1' => {"vanilla" => 0.75, "chocolate" => 0.25},
       	'Bowl 2' => {"vanilla" => 0.5, "chocolate" => 0.5}
	}

	def likelihood(data, hypo)
        mix = @@mixes[hypo]
        like = mix[data]
        return like
	end

	def update(data)
        @d.each_key{|hypo|
            like = self.likelihood(data, hypo)
            self.mult(hypo, like)
        }
        self.normalize()
	end

end

hypos = ['Bowl 1', 'Bowl 2']
pmf = Cookie.new(hypos)
pmf.update('vanilla')
pmf.update('vanilla')
pmf.update('chocolate')

pmf.print_dict
p pmf
