#################################################
##    Bayes Utility
##    2014-11-
##    Python版書籍からの移植版
#################################################
require 'numo/narray'
require 'numo/gnuplot'
require 'numo/narray'
require 'randomext'
require './kde'

class Interpolator
    # Represents a mapping between sorted sequences; performs linear interp.
    # 線形補間クラス、点と点との中間の点を計算する
	#   xs,ys: sorted list
    def initialize(xs, ys)
        @xs = xs
        @ys = ys
	end

    def lookup(x)
        # Looks up x and returns the corresponding value of y.
        return bisect(x, @xs, @ys)
	end

    def reverse(y)
        # Looks up y and returns the corresponding value of x.
        return bisect(y, @ys, @xs)
	end

    def bisect(x, xs, ys)
        return ys[0]  if x <= xs[0]
        return ys[-1] if x >= xs[-1]
        i = xs.bsearch_right(x)
        frac = 1.0 * (x - xs[i - 1]) / (xs[i] - xs[i - 1])
        y = ys[i - 1] + frac * 1.0 * (ys[i] - ys[i - 1])
        return y
	end
end

class PmfSuite
	# Represents a probability mass function.
	# Values can be any hashable type; probabilities are floating-point.
	# Pmfs are not necessarily normalized.
	include Enumerable
	attr_accessor :name 
    def initialize(values = nil, name: nil)
		@hypos			# sequence of hypotheses
        @d = {}
		@log = false	# flag whether the distribution is under a log transform
        @name = name
        
        return   if values.nil?

        init_methods = [
            method(:init_pmf),
            method(:init_hash),
            method(:init_sequence),
            method(:init_failure)
            ]

        init_methods.each{|f|
            begin
                f.call(values)
                break
            rescue
            	next
            end
		}
	end

    def init_sequence(values)
    	# Initializes with a sequence of equally-likely values.
		raise 'data unmatch' unless values.is_a?(Array) || values.is_a?(String)
		values = values.split("") if values.is_a?(String)
        values.each{|v| incr(v)}
	end

    def init_pmf(values)
    	# Initializes with a Pmf.
		raise 'data unmatch' unless values.is_a?(PmfSuite)
        values.get_dict().each{|value, prob| set(value, prob)}
	end
	
    def init_hash(values)
   		# Initializes with a hash.
		raise 'data unmatch' unless values.is_a?(Hash)
        values.each{|value, prob| set(value, prob)}
	end

    def init_failure(values)
        @log = 'None of the initialization methods worked.'
        raise 'None of the initialization methods worked.'
	end

    def len()
        return @d.size
	end

	def each
    	@d.each do |x|
      	yield x
    	end
  	end

    def copy
    	# Returns a deep copy of PmfSuite.  深いコピー
        new_d = Marshal.load(Marshal.dump(self))   
        return new_d
	end

    def scale(factor)
    	# Multiplies the values by a factor.
        new_d = @d.each{|val, prob|
            [val * factor, prob]
        }.to_h
        return PmfSuite.new(new_d)
	end
	
    def log(m = nil)
    	# Log transforms the probabilities.
        # Removes values with probability 0.
		# Normalizes so that the largest logprob is 0.
        raise "PmfSuite already under a log transform"  if @log_log
        @log_log = true
        m = max_like()  if m.nil?
        self.each{|x, prob|
            if prob
                set(x, Math.log(prob / m))
            else
                remove(x)
            end
        }
        return self
	end		
			
    def exp(m = nil)
   		# Exponentiates the probabilities.
		#   m: how much to shift the ps before exponentiating
		# If m is None, normalizes so that the largest prob is 1.
        m = max_like()  if m.nil?
        each{|x, y| set(x, Math.exp(y - m))}
        return self
	end

    def get_dict()
    	# オリジナルの@dと同じオブジェクトを返す。中身のコピーではない。
        return @d.clone
	end

    def set_dict(d)
        @d = d
	end

    def values()
    	# Gets an unsorted sequence of values.
		#   Note: 頻度／確率値ではなくPmfの値(説明変数)リストを取得
		return @d.keys()
	end

    def items()
        return @d.to_a
	end

    def print_dict
        @d.each{|val, prob|  print val," = ", prob, "\n"}
	end

    def set(x, y = 0)
    	# Sets the freq/prob associated with the value x.
		#   x: number value
    	#   y: number freq or prob
        @d[x] = y
	end

    def incr(x, term = 1)
    	# Increments the freq/prob associated with the value x.
		#   x: number value
    	#   term: how much to increment by
        @d[x] ? @d[x] += term : @d[x] = term
    end

    def mult(x, factor)
    	# Scales the freq/prob associated with the value x.
		#   x: number value
    	#   factor: how much to multiply by
       	return @d[x] *= factor	if @d[x]
		return @d[x.to_i] *= factor	if @d[x.to_i]
		return @d[x.to_f] *= factor	if @d[x.to_f]
	end

    def remove(x)
        @d.delete(x)
	end

    def total
        return @d.values.inject(:+)
	end

    def max_like
    	# Returns the largest frequency/probability in the map.
        return (@d.max{ |x, y| x[1] <=> y[1] })[1]
	end

#####  以降確率操作関連メソッド、かつてのPmfクラス
    def prob(x)
    	# Gets the probability associated with the value x.
		#   x: number value
		#   Returns: float probability
		return @d[x]  		if @d[x]
		return @d[x.to_i]	if @d[x.to_i]
		return @d[x.to_f]	if @d[x.to_f]
        return 0	# default
	end

    def probs(xs)
        return xs.collect{|x| prob(x)}
	end

    def make_cdf
        return make_cdf_from_pmf(self)
	end

    def prob_greater(x)
    	t = items.select{|val, plob| val > x}
        return t.transpose()[1].sum()
	end

    def prob_less(x)
        t = items.select{|val, plob| val < x}
        return t.transpose()[1].sum()
	end

    def normalize(fraction = 1.0)
    	# Normalizes this PMF so the sum of all probs is fraction.
		#   fraction: what the total should be after normalization
		#   Returns: the total probability before normalizing
		raise "Pmf is under a log transform// #{@log}" if @log
        total = total()
        if total == 0
        	@log = "Normalize: total probability is zero."
        	raise "total probability is zero."
        end
        @d.each{ |x, y| @d[x] = y.to_f * fraction.to_f / total }
        return total
	end

    def random()
    	# Chooses a random element from this PMF.
		#   Returns: float value from the Pmf
		raise 'Pmf contains no values.'   if len == 0
		target = rand(total())
        total = 0.0
        @d.each{|val, prob|
			total += prob
			return val  if total > target
		}
		return false
	end

    def mean()
    	# Computes the mean of a PMF.
		#   Returns: float mean
        mu = 0.0
        @d.each{|hypo, prob| mu += hypo * prob}
        return mu
	end

    def var(mu = nil)
    	# Computes the variance of a PMF.
		#   mu: the point around which the variance is computed;
		#		 if omitted, computes the mean
		#   Returns: float variance
        mu = mean()  if mu == nil
        var = 0.0
        @d.each{ |x, y| var += y * (x - mu) ** 2 }
        return var
	end

    def maximum_likelihood()
    	# Returns the value with the highest probability.
		#   Returns: float probability
         return (@d.max{ |x, y| x[1] <=> y[1] })[0]
	end

    def credible_interval(percentage = 90)
    	# Computes the central credible interval.
		#    Returns: sequence of two floats, low and high
        cdf = make_cdf()
        return cdf.credible_interval(percentage)
 	end

    def +(other)
    	# Computes the Pmf of the sum of values drawn from self and other.
        begin
            return add_pmf(other)
        rescue
            return add_constant(other)
        end
 	end

    def add_pmf(other)
    	# Computes the Pmf of the sum of values drawn from self and other.
        pmf = PmfSuite.new()
        items.each{|v1, p1|
        	other.items.each{|v2, p2|
        		pmf.incr(v1 + v2, p1 * p2)
        	}
        }
        pmf.normalize
        return pmf
	end

    def add_constant(other)
    	# Computes the Pmf of the sum a constant and  values from self.
        pmf = PmfSuite.new()
        items.each{|v1, p1|
        	pmf.set(v1 + other, p1)
        }
#        pnf.normalize
        return pmf
	end

    def -(other)
    	# Computes the Pmf of the diff of values drawn from self and other.
        pmf = PmfSuite.new()
        items.each{|v1, p1|
        	other.items.each{|v2, p2|
        		pmf.incr(v1 - v2, p1 * p2)
        	}
        }
        pmf.normalize
        return pmf
 	end
	
	def percentile(percentage = 90)
		# Computes a percentile of a given Pmf.
		#    percentage: float 0-100
    	x = percentage / 100.0
    	total = 0.0
    	items.each{|val, prob|
    		total += prob
        	return val  if total >= x
    	}
	end

    def max(k)
    	# Computes the CDF of the maximum of k selections from this dist.
    	# 個人的には意味不明
		#   k: int
		#   returns: new Cdf
        cdf = make_cdf
        return cdf.max(k)
	end

	def render
		# plot formatの生成
		dd = items.transpose
		return [dd[0], dd[1], {w:'lines', t: @name}]		
	end

	def make_plot(x_label: "x_label", y_label: "PMF", title: "title")
		make_plot_2d(render, x_label: x_label, y_label: y_label, title: title, grid: "on", )
	end

# Represents a suite of hypotheses and their probabilities.
# Suite用メソッド、かつてのSuiteクラス
    def update(data)
    	# Updates each hypothesis based on the data.
		#   data: any representation of the data
		#   returns: the normalizing constant
		@d.each{|k, v| @d[k] = mult(k, likelihood(data, k)) }
        return normalize
	end

    def log_update(data)
    	# Updates a suite of hypotheses based on new data.
		# Modifies the suite directly; if you want to keep the original, make a copy.
		# Note: unlike Update, LogUpdate does not normalize.
		#   data: any representation of the data
        values().each{|hypo|
            like = log_likelihood(data, hypo)
            incr(hypo, like)
        }
	end

    def update_set(dataset)
    	# Updates each hypothesis based on the dataset.
		# This is more efficient than calling Update repeatedly because
		# it waits until the end to Normalize.
		# Modifies the suite directly; if you want to keep the original, make a copy.
		#   dataset: a sequence of data
		#   returns: the normalizing constant
        dataset.each{|data|
            values.each{|hypo|
                like = likelihood(data, hypo)
                mult(hypo, like)
            }
        }
        normalize
        return self
	end


    def log_update_set(dataset)
        # Updates each hypothesis based on the dataset.
        # Modifies the suite directly; if you want to keep the original, make
        # a copy.
        #   dataset: a sequence of data
        #   returns: None
        dataset.each{|data| log_update(data) }
	end

	def likelihood
			raise "no method" #"liklihoodメソッドの未定義、定義して下さい。\n"
	end
	def log_likelihood
			raise "no method" #"liklihoodメソッドの未定義、定義して下さい。\n"
	end

    def probs_to_odds()		# make_odds
		# Values with prob=0 are removed.
        self.items.each{|hypo, prob|
        	if prob
                self.ret(hypo, odds(prob))
            else
                self.remove(hypo)
			end
		}
	end

    def odds_to_probs()		# make_probs
        self.items.each{|hypo, odds| self.set(hypo, probability(odds)) }
	end
end

class Joint < PmfSuite
    # 数字組配列 と 頻度/確率のセットを表すクラス　(joint distribution)
    def marginal(i)
        # Gets the marginal distribution of the indicated variable.
        # 変数 i についての周辺分布を作る
       pmf = PmfSuite.new()
        items.each{|vs, prob|
            pmf.incr(vs[i], prob)
        }
        return pmf
	end

    def conditional(i, j, val)
        # Gets the conditional distribution of the indicated variable.
        # 条件付き分布を作る
		# Distribution of vs[i], conditioned on vs[j] == val.
		#     i: index of the variable we want
        #     j: which variable is conditioned on
        #   val: the value the jth variable has to have
		#   Returns: Pmf

        pmf = PmfSuite.new(name: val.to_s)
        items.each{|vs, prob|
            next if vs[j] != val
            pmf.incr(vs[i], prob)
		}
        pmf.normalize()
        return pmf
	end

    def max_like_interval(percentage = 90)
		# percentage=90の場合は、高い方からの累積確率がTOTAL90% までに入る要素集合を返す
		#     Returns: list of values from the suite
        interval = []
        total = 0

        t = items().to_a 
        t.sort!{|a, b| a[1] <=> b[1]}.reverse!()

        t.each{|val, prob|
            interval << val
            total += prob
            break if total >= percentage / 100.0
		}
        return interval
	end
	
	def render
		# plot formatの生成
		dd = items.each.collect{|key, value| [key[0], key[1], value] }.transpose
		return [dd[0], dd[1], dd[2], {w:'pm3d at s', t: @name}]		
	end

	def make_plot(x_label: "x_label", y_label: "y_label", title: "title")
		make_plot_3d(render, x_label: x_label, y_label: y_label, title: title, contour: "both")
	end
end

def make_joint(pmf1, pmf2)
    # Joint distribution of values from pmf1 and pmf2.
	# Returns: Joint pmf of value pairs
    joint = Joint.new()
    pmf1.items.each {|v1, p1|
        pmf2.items.each{|v2, p2|
        	joint.set([v1, v2], p1 * p2)
        }
    }
    return joint
end

def make_suite_from_list(t, name: nil)
	#   t: unsorted sequence of numbers
    suite = PmfSuite.new(t, name: name)
    suite.normalize()
    return suite
end

def make_suite_from_dict(d, name: nil)
    suite = PmfSuite.new(d, name: name)
    suite.normalize()
    return suite
end

def make_suite_from_cdf(cdf, name: nil)
    # Makes a normalized Suite from a Cdf object.
   suite = PmfSuite.new(name: name)
    prev = 0.0
    cdf.items.each{|val, prob|
        suite.incr(val, prob - prev)
        prev = prob
	}
    return suite
end

def make_suite_from_items(t, name: nil)
	# Makes a PMF from a sequence of value-probability pairs
	suite = PmfSuite.new(t.to_h, name: name)
    suite.normalize()
    return suite
end

def make_mixture(metapmf, name: nil)
	# Make a mixture distribution.
	#   metapmf: Pmf that maps from Pmfs to probs.
    mix = PmfSuite.new()
    metapmf.each{|pmf, p1|
    	pmf.each{|x, p2|
    		mix.incr(x, p1 * p2)
    	}
    }
    return mix
end

def pmf_prob_less(pmf1, pmf2)
    # Probability that a value from pmf1 is less than a value from pmf2.
	# Returns: float probability　　全組み合わせの中でpmf1の方が小さい確率
    total = 0.0
    pmf1.items.each{|v1, p1| pmf2.items.each{|v2, p2| total += p1 * p2  if v1 < v2 }}
    return total
end

def pmf_prob_greater(pmf1, pmf2)
    # Probability that a value from pmf1 is less than a value from pmf2.
	# Returns: float probability　　全組み合わせの中でpmf1の方が小大きい確率
    total = 0.0
    pmf1.items.each{| v1, p1 | pmf2.items.each{| v2, p2 | total += p1 * p2  if v1 > v2 }}
    return total
end

def pmf_prob_equal(pmf1, pmf2)
    # Probability that a value from pmf1 equals a value from pmf2.
    # Returns: float probability　　全組み合わせの中で、pmf1/2が同じになる確率
    total = 0.0
    pmf1.items.each{| v1, p1 | pmf2.items.each{| v2, p2 | total += p1 * p2 if v1 == v2 }}
    return total
end

def make_uniform_pmf(low, high, n)
	# Make a uniform Pmf.
	#   low: lowest value (inclusive)
	#   high: highest value (inclusize)
	#   n: number of values
    pmf = PmfSuite.new()
    n = (high - low) / (n - 1)
    low.step(high, n){|x| pmf.set(x, 1)}
    pmf.normalize()
    return pmf
end

def random_sum(dists)
	# Chooses a random value from each dist and returns the sum.
	#   dists: sequence of Pmf or Cdf objects
	#   returns: numerical sum
    return dists.collect{|dist| dist.random()}.sum
end

def sample_sum(dists, n)
	# Draws a sample of sums from a list of distributions.
	#   dists: sequence of Pmf or Cdf objects
	#   n: sample size
	#   returns: new Pmf of sums
    return make_suite_from_list(n.times.collect{random_sum(dists)})
end

def random_max(dists)
	return dists.collect{|dist| dist.random}.max
end

def sample_max(dists, n)
    return make_suite_from_list(n.times.collect{random_max(dists)})
end

class Cdf
	# Represents a cumulative distribution function.
	#     @xs: sequence of values
	#     @ps: sequence of probabilities
	attr_reader  :xs, :ps
	attr_accessor	:name
    def initialize(xs = nil, ps = nil, name: nil)
    	raise "Cdf data not Array" unless xs.is_a?(Array) && ps.is_a?(Array)
        @xs = (xs.nil? ?  [] : xs)
        @ps = (ps.nil? ?  [] : ps)
        @name = name
	end

	def copy(name: nil)
		# 浅いコピー、実体は同じもの
        return Cdf.new(@xs.clone, @ps.clone, name: name)
	end

    def make_pmf(name: nil)
        return make_suite_from_cdf(self, name: name)
	end

    def values()
        return @xs.sort
	end

    def items()
    	return @xs.zip(@ps)
	end

    def append(x, y)
    	# Add an (x, y) pair to the end of this CDF.
		# 既存のCDFへの追加は推奨しない。ゼロからCDFを作る場合に使用。正しいCDFになっているかは未保証
		@xs.push(x)
		@ps.push(y)
	end

    def shift(term)
        new_xs = @xs.map{|x| x + term}
        return Cdf.new(new_xs, @ps)
	end

    def scale(factor)
    	# Multiplies the xs by a factor.
        new_xs = @xs.map{|x| x * factor.to_f}
        return Cdf.new(new_xs, @ps)
	end

    def prob(x)
     	# Returns CDF(x), xの値に対応する累積確率を返す
		return 0.0 if x < @xs[0]
		index = @xs.bsearch_right(x)
        return @ps[index - 1]
	end

    def value(x)
    	# Returns InverseCDF(p), the value that corresponds to probability p.
		# 	p: number in the range [0, 1]
		#   Returns: number value
        raise 'Probability p must be in range [0, 1]'   if x < 0 or x > 1
        return @xs[0]	if x == 0 
        return @xs.last	if x == 1
		index = @ps.bsearch_right(x)
        if x == @ps[index - 1]
            return @xs[index - 1]
        else
            return @xs[index]
		end
	end

    def percentile(x)
    	# Returns the value that corresponds to percentile p.
        return value(x / 100.0)
	end

    def random
    	# Chooses a random value from this distribution.
        return value(rand)
	end

    def sample(n)
    	# Generates a random sample from this distribution.
		ar = []
		n.to_i.times{ ar << random} 
        return ar
	end

    def mean()
    	# Computes the mean of a CDF.
        old_p = 0
        total = 0.0
        @xs.zip(@ps).each{|x, new_p|
            tmp = new_p - old_p
            total += tmp * x
            old_p = new_p
        }
        return total
	end

    def credible_interval(percentage = 90)
    	# Computes the central credible interval.
		#   Returns: sequence of two floats, low and high
        prob = (1 - percentage / 100.0) / 2
        interval = value(prob), value(1 - prob)
        return interval
	end

    def max(k)
    	# Computes the CDF of the maximum of k selections from this dist.
		#   k: int
		#   returns: new Cdf
        new_ps = @ps.map{|x| x**k}
        cdf = Cdf.new(@xs.clone, new_ps)
        return cdf
	end

	def render
		# plot formatの生成
		return [@xs, @ps, {w:'lines', t: @name}]
=begin
    if transform == 'exponential':
        complement = True
        scale['yscale'] = 'log'

    if transform == 'pareto':
        complement = True
        scale['yscale'] = 'log'
        scale['xscale'] = 'log'

    if complement:
        ps = [1.0-p for p in ps]

    if transform == 'weibull':
        xs.pop()
        ps.pop()
        ps = [-math.log(1.0-p) for p in ps]
        scale['xscale'] = 'log'
        scale['yscale'] = 'log'

    if transform == 'gumbel':
        xs.pop(0)
        ps.pop(0)
        ps = [-math.log(p) for p in ps]
        scale['yscale'] = 'log'
=end
	end
	
	def make_plot(x_label: "x_label", y_label: "CDF", title: "title")
		make_plot_2d(render, x_label: x_label, y_label: y_label, title: title, grid: "on", )
	end
end

def make_cdf_from_items(items, name: nil)
	# Makes a cdf from an unsorted sequence of (value, frequency) pairs.
	#	items: unsorted sequence of (value, frequency) pairs
	#   name: string name for this CDF
	#   Returns: cdf -- list of (value, fraction) pairs
    runsum = 0
    xs = []
    cs = []

    items.sort.each{|value, count|
        runsum += count.to_f
        xs << value.to_f
        cs << runsum
	}
    total = runsum.to_f
    ps = cs.map{|c| c / total}
    cdf = Cdf.new(xs, ps)
    return cdf
end

def make_cdf_from_dict(d, name: nil)
	# Makes a CDF from a dictionary that maps values to frequencies.
    return make_cdf_from_items(d.to_a)
end

def make_cdf_from_pmf(pmf, name: nil)
	# Makes a CDF from a Pmf object.
   return make_cdf_from_items(pmf.items, name: pmf.name)
end

def make_cdf_from_list(seq, name: nil)
	# Creates a CDF from an unsorted sequence.
    pmf = PmfSuite.new(seq.collect{|x| x.to_f}.sort)
    return make_cdf_from_pmf(pmf)
end

class Pdf
    def density(x)
		#  Returns: float probability density
		raise "UnimplementedMethod"
	end

    def make_pmf(xs)
        pmf = PmfSuite.new()
        xs.each{|x| pmf.set(x, density(x)) }
        pmf.normalize()
        return pmf
	end
end

class GaussianPdf < Pdf
	def initialize(mu, sigma, name: nil)
		#  mu: mean
		#  sigma: standard deviation
        @mu = mu
        @sigma = sigma
        @name = name
	end

    def density(x)
        return eval_gaussian_pdf(x, @mu, @sigma)
	end
end

class EstimatedPdf < Pdf
     def initialize(sample, bw_width = "silverman", name: nil)
        @kde = GaussianKde.new(sample, bw_width)
        @name = name
	end

	def resample(num)
		@kde.resample(num)
	end

    def density(x)
        return @kde.evaluate(x)
	end

    def make_pmf(xs)
        ps = @kde.evaluate(xs)
        pmf = make_suite_from_items(xs.zip(ps))
        return pmf
	end
end

def eval_gaussian_pdf(x, mu, sigma)
	# ガウス分布からの確率値を計算
    return Math::E ** (-(x - mu)**2/2/(sigma**2)) / (Math.sqrt(2 * Math::PI) * sigma)
end

def make_gaussian_pmf(mu, sigma, num_sigmas, n = 201)
	# num_sigmas: how many sigmas to extend in each direction
    pmf = PmfSuite.new()
    low = mu.to_f - num_sigmas.to_f * sigma.to_f
    high = mu.to_f + num_sigmas.to_f * sigma.to_f
	n = (high - low) / (n - 1)

    low.step(high, n){ |x|
        p = eval_gaussian_pdf(x, mu, sigma)
        pmf.set(x, p)
    }
    pmf.normalize
    return pmf
end

def eval_poisson_pmf(k, lam)
	# ポアソンPMFからの確率値を計算
	# k: number of events
	# lam: parameter lambda in events per unit time
	factr(k)
    return lam.to_f ** k * Math.exp(-lam.to_f) / factr(k)
end

def make_poisson_pmf(lam, high, st = 1)
	# lam: parameter lambda in events per unit time
	# high: upper bound of the Pmf
    pmf = PmfSuite.new()
    sigma = Math.sqrt(lam)
	0.step(high + 1, st) { |k|
		if lam > 20
			pro = eval_gaussian_pdf(k, lam, sigma)
		else
        	pro = eval_poisson_pmf(k, lam)
        	pro = 0 if pro.nan?
		end
        pmf.set(k, pro)
    }
    pmf.normalize
    return pmf
end

def eval_exponential_pdf(x, lam)
	# 指数分布からの確率値を計算
	# x: value
	# lam: parameter lambda in events per unit time
	return lam * Math.exp(-lam * x)
end

def eval_exponential_cdf(x, lam)
	# Evaluates CDF of the exponential distribution with parameter lam."""
	return 1 - Math.exp(-lam * x)
end

def make_exponential_pmf(lam, high, n=201)
	# lam: parameter lambda in events per unit time
	# high: upper bound
    pmf = PmfSuite.new()
    n = (high.to_f - 0) / n
    0.step(high, n) { |x|
        p = eval_exponential_pdf(x, lam)
        pmf.set(x, p)
    }
    pmf.normalize
    return pmf
end

def standard_gaussian_cdf(x, root2 = Math.sqrt(2))
	# Evaluates the CDF of the standard Gaussian distribution.
	#    See http://en.wikipedia.org/wiki/Normal_distribution
	#    #Cumulative_distribution_function
	# x: float
	# Returns: float
    return (Math.erf(x / root2) + 1) / 2
end

def gaussian_cdf(x, mu = 0, sigma = 1)
	# Evaluates the CDF of the gaussian distribution.
	# x: float
	# mu: mean parameter
	# sigma: standard deviation parameter
	# Returns: float
    return standard_gaussian_cdf((x - mu).to_f / sigma)
end

=begin
def gaussian_cdf_inverse(p, mu = 0, sigma = 1)
    # Evaluates the inverse CDF of the gaussian distribution.
	#  See http://en.wikipedia.org/wiki/Normal_distribution#Quantile_function  
	#   p: float
	#   mu: mean parameter
    #   sigma: standard deviation parameter
    #   Returns: float
    x = root2 * erfinv(2 * p - 1)
    return mu + x * sigma
end
=end


class Beta
    # Represents a Beta distribution.
	# See http://en.wikipedia.org/wiki/Beta_distribution
	attr_accessor :name
	
    def initialize(alpha = 1, beta = 1, name: nil)
        @alpha = alpha
        @beta = beta
        @name = name
	end

    def update(data)
		# data: pair of int (heads, tails)
        heads, tails = data
        @alpha += heads
        @beta += tails
	end

    def mean
        return @alpha.to_f / (@alpha + @beta)
	end

    def random
        return random.betavariate(@alpha, @beta)
	end

    def eval_pdf(x)
        return x ** (@alpha - 1) * (1 - x) ** (@beta - 1)
	end

    def make_pmf(steps = 101)
		# alpha/betaが1未満の場合は両端(x=0or1)が無限大になるので注意必要
        # その場合は一旦CDFに変換して評価する
        if @alpha < 1 or @beta < 1
            cdf = make_cdf()
            pmf = cdf.make_pmf()
            return pmf
        end

        xs = 0.step(1.0, 1 / (steps - 1.0)).to_a
        probs = xs.collect{|x| eval_pdf(x)}
        pmf = make_suite_from_dict(xs.zip(probs).to_h, name: @name)
        return pmf
	end

    def make_cdf(steps = 101)
        xs = 0.step(1.0, 1 / (steps - 1.0)).to_a
        ps = xs.collect{|x| p_beta(x, @alpha, @beta)}
        cdf = Cdf.new(xs, ps)
        return cdf
	end
end

class Dirichlet
	#   See https://ja.wikipedia.org/wiki/%E3%83%87%E3%82%A3%E3%83%AA%E3%82%AF%E3%83%AC%E5%88%86%E5%B8%83
	attr_reader  :n, :params
    def initialize(n, conc = 1.0, name:nil)
        #   n: number of dimensions
        #   conc: concentration parameter (smaller yields more concentration)
        raise 'A Dirichlet distribution with n<2 makes no sense'  if n < 2
        @n = n
        @params = Numo::DFloat.ones(n) * conc
        @name = name
	end

    def update(data)
		#   data: sequence of observations, in order corresponding to params
        @params[0..data.size - 1] += data
	end

    def random
		#   Returns: normalized vector of fractions
		a = Numo::DFloat.new(@n)
		a.store( @params.to_a.collect{|a| Random::DEFAULT.gamma(a)} )
		return a / a.sum
	end

    def likelihood(data)
        # Computes the likelihood of the data.
		# Selects a random vector of probabilities from this distribution.
		# Returns: float probability
        return 0  if @n < data.size
        pr = self.random()
        q = pr[0..data.size-1] ** data
        return q.prod
	end

    def log_likelihood(data)
        # Computes the log likelihood of the data.
		# Selects a random vector of probabilities from this distribution.
		# Returns: float log probability
        return - Float::INFINITY  if @n < data.len

        x = self.random()
        y = Numo::NMath.log(x[0..data.size-1]) * data
        return y.sum
	end

    def marginal_beta(i)
        # i番目要素についての周辺分布を計算する
		# See http://en.wikipedia.org/wiki/Dirichlet_distribution
		#   Returns: Beta object
        alpha0 = @params.sum()
        alpha = @params[i]
        return Beta.new(alpha, alpha0 - alpha)
	end
	
    def predictive_pmf(xs, name: nil)
        # 予測分布の作成
		#   xs: values to go into the Pmf
		#   Returns: Pmf that maps from x to the mean prevalence of x
        alpha0 = @params.sum()
        ps = @params / alpha0
        return make_suite_from_items(xs.zip(ps.to_a), name: name)
	end
end

def jitter(values, jitter = 0.5)
    # Jitters the values by adding a uniform variate in (-jitter, jitter).
    return values.collect{|x| x + rand(-jitter..jitter)}
end

#########################################################################
###
###　  　　 以下　Array クラスの拡張
###
#########################################################################
class Array
  def sum
    self.reduce{|sum, n| sum.to_f + n.to_f}
  end

  def mean
    self.sum / self.size
  end

  def var
    m = self.mean
    self.reduce(0){|sum,b| sum.to_f + (b.to_f - m) ** 2 } / (self.size - 1)
  end

  def sd
    Math.sqrt(self.var)
  end

####### 以下　python bisect のruby実装。　ソート済み数列に大小関係を崩さずに挿入できる挿入点を探す
  def insert_right(item, lo=0, hi=nil)
    if lo < 0
        raise IndexError.new('lo must be non-negative')
    end
    if hi == nil
        hi = self.length
    end
    while lo < hi
        mid = ((lo+hi)/2).truncate
        if item < self[mid]
          hi = mid
        else 
          lo = mid+1
        end
    end
    self.insert(lo, item)
  end

  def bsearch_right(item, lo=0, hi=nil)
    if lo < 0
        raise IndexError.new('lo must be non-negative')
    end
    if hi == nil
        hi = self.length
    end
    while lo < hi
        mid = ((lo+hi)/2).truncate
        if item < self[mid]
          hi = mid
        else
          lo = mid+1
        end
    end
    return lo
  end
  
  def insert_left(item, lo=0, hi=nil)
    if lo < 0
        raise IndexError.new('lo must be non-negative')
    end
    if hi == nil
        hi = self.length
    end
    while lo < hi
        mid = ((lo+hi)/2).truncate
        if self[mid] < item 
          lo = mid+1
        else
          hi = mid
        end
    end
    self.insert(lo, item)
  end

  def bsearch_left(item, lo=0, hi=nil)
    if lo < 0
        raise IndexError.new('lo must be non-negative')
    end
    if hi == nil
        hi = self.length
    end
    while lo < hi
        mid = ((lo+hi)/2).truncate
        if self[mid] < item
          lo = mid+1
        else
          hi = mid
        end
    end
    return lo
  end

################### 以下一次の最小二乗近似式の計算　########################
  def reg_line(y)
    # 以下の場合は例外スロー
    # - 引数の配列が Array クラスでない
    # - 自身配列が空
    # - 配列サイズが異なれば例外
    raise "Argument is not a Array class!"  unless y.class == Array
    raise "Self array is nil!"              if self.size == 0
    raise "Argument array size is invalid!" unless self.size == y.size

    # x^2 の総和
    sum_xx = self.inject(0) { |s, a| s += a * a }

    # x * y の総和
    sum_xy = self.zip(y).inject(0) { |s, a| s += a[0] * a[1] }

    # 切片 a
    a  = sum_xx * y.sum - sum_xy * self.sum
    a /= (self.size * sum_xx - self.sum * self.sum).to_f

    # 傾き b
    b  = self.size * sum_xy - self.sum * y.sum
    b /= (self.size * sum_xx - self.sum * self.sum).to_f

    return {intercept: a, slope: b}
  end

########################### 以下相関係数の計算 ##############################
  def corr(y)
    # 以下の場合は例外スロー
    # - 引数の配列が Array クラスでない
    # - 自身配列が空
    # - 配列サイズが異なれば例外
    raise "Argument is not a Array class!"  unless y.class == Array
    raise "Self array is nil!"              if self.size == 0
    raise "Argument array size is invalid!" unless self.size == y.size

    # x の相加平均, y の相加平均 (arithmetic mean)
    mean_x = self.mean
    mean_y = y.mean

    # x と y の共分散 (covariance)
    cov = self.zip(y).inject(0) { |s, a| s += (a[0] - mean_x) * (a[1] - mean_y) } / self.size

    # 相関係数 (correlation coefficient)
    return cov / self.sd / y.sd
  end
##########################  系列相関係数の計算　##############################
  def serial_corr
    # 自身配列が空の場合は例外
    raise "Self array is nil!"  if self.size == 0

	# 時系列の一次近似式の計算
	x = 1.upto(self.size).to_a
	ans = x.reg_line(self)
	
	# 残差の計算　→Array
	e = x.each_index.collect{|i| self[i] - (ans[:intercept] + x[i] * ans[:slope])}

    # 1データズレの場合積和計算
    sumproduct = 0
    e[1..-1].each_index{|i| sumproduct += e[i] * e[i - 1]}

    # x と y の共分散の分子 (covariance)
    sumsq = e[1..-1].inject{|sum, a| sum += a**2}

    # 系列相関係数 (Serial correlation coefficient)
    return sumproduct / sumsq
  end
end
=begin
  ****  以下　Pyhtonルーチンをそのままコピー移植未完成  ****
def cov(xs, ys, mux = nil, muy = nil)
    # Computes Cov(X, Y).
	#   xs: sequence of values
    #   ys: sequence of values
    #   mux: optional float mean of xs
    #   muy: optional float mean of ys
	#	Returns:　Cov(X, Y)
    mux = thinkstats.Mean(xs)  if mux.nil?
    muy = thinkstats.Mean(ys)  if muy.nil?

    total = 0.0
    for x, y in zip(xs, ys):
        total += (x-mux) * (y-muy)

    return total / len(xs)
end

def Corr(xs, ys):
    """Computes Corr(X, Y).

    Args:
        xs: sequence of values
        ys: sequence of values

    Returns:
        Corr(X, Y)
    """
    xbar, varx = thinkstats.MeanVar(xs)
    ybar, vary = thinkstats.MeanVar(ys)

    corr = Cov(xs, ys, xbar, ybar) / math.sqrt(varx * vary)

    return corr


def SerialCorr(xs):
    """Computes the serial correlation of a sequence."""
    return Corr(xs[:-1], xs[1:])


def SpearmanCorr(xs, ys):
    """Computes Spearman's rank correlation.

    Args:
        xs: sequence of values
        ys: sequence of values

    Returns:
        float Spearman's correlation
    """
    xranks = MapToRanks(xs)
    yranks = MapToRanks(ys)
    return Corr(xranks, yranks)


def LeastSquares(xs, ys):
    """Computes a linear least squares fit for ys as a function of xs.

    Args:
        xs: sequence of values
        ys: sequence of values

    Returns:
        tuple of (intercept, slope)
    """
    xbar, varx = thinkstats.MeanVar(xs)
    ybar, vary = thinkstats.MeanVar(ys)

    slope = Cov(xs, ys, xbar, ybar) / varx
    inter = ybar - slope * xbar

    return inter, slope


def FitLine(xs, inter, slope):
    """Returns the fitted line for the range of xs.

    xs: x values used for the fit
    slope: estimated slope
    inter: estimated intercept
    """
    fxs = min(xs), max(xs)
    fys = [x * slope + inter for x in fxs]
    return fxs, fys


def Residuals(xs, ys, inter, slope):
    """Computes residuals for a linear fit with parameters inter and slope.

    Args:
        xs: independent variable
        ys: dependent variable
        inter: float intercept
        slope: float slope

    Returns:
        list of residuals
    """
    res = [y - inter - slope*x for x, y in zip(xs, ys)]
    return res


def CoefDetermination(ys, res):
    """Computes the coefficient of determination (R^2) for given residuals.

    Args:
        ys: dependent variable
        res: residuals
        
    Returns:
        float coefficient of determination
    """
    ybar, vary = thinkstats.MeanVar(ys)
    resbar, varres = thinkstats.MeanVar(res)
    return 1 - varres / vary


def MapToRanks(t):
    """Returns a list of ranks corresponding to the elements in t.

    Args:
        t: sequence of numbers
    
    Returns:
        list of integer ranks, starting at 1
    """
    # pair up each value with its index
    pairs = enumerate(t)
    
    # sort by value
    sorted_pairs = sorted(pairs, key=lambda pair: pair[1])

    # pair up each pair with its rank
    ranked = enumerate(sorted_pairs)

    # sort by index
    resorted = sorted(ranked, key=lambda trip: trip[1][0])

    # extract the ranks
    ranks = [trip[0]+1 for trip in resorted]
    return ranks


def CorrelatedGenerator(rho):
    """Generates standard normal variates with correlation.

    rho: target coefficient of correlation

    Returns: iterable
    """
    x = random.gauss(0, 1)
    yield x

    sigma = math.sqrt(1 - rho**2);    
    while True:
        x = random.gauss(x * rho, sigma)
        yield x


def CorrelatedNormalGenerator(mu, sigma, rho):
    """Generates normal variates with correlation.

    mu: mean of variate
    sigma: standard deviation of variate
    rho: target coefficient of correlation

    Returns: iterable
    """
    for x in CorrelatedGenerator(rho):
        yield x * sigma + mu

=end
#########################################################################
###
### String クラスの拡張  文字列が数字か確認する
###
#########################################################################
class String
  def int_valid?
    Integer(self)
    true
  rescue ArgumentError
    false
  end

  def float_valid?
    Float(self)
    true
  rescue ArgumentError
    false
  end
end

###################################################################
###    Randomクラスの拡張　
###    乱数の拡張は　Randomext gemで対応。以下の拡張ができる
###		normal (Gaussian), lognormal, Cauthy, levy, exponential, Laplace, Rayleigh
###		Weibull, Gumbel, gamma, beta, power, Chi-square, F, t, Wald (inverse Gaussian)
###		Pareto, logistic, von Mises, Non-Central Chi-Square, Non-Central t, Planck
###		Bernoulli, binomial, Poisson, geometric, negative binomial, log series
###		Zipf-Mandelbrot, zeta
###　　　アルゴリズムの出典は
###  	Almost all algorithms are based on: 四辻哲章, "計算機シミュレーションのための確率分布乱数生成法", プレアデス出版 (2010)
###   使用例は
###		require 'randomext'
###		random_numbers = Array.new(100){ Random::DEFAULT.normal(0.0, 2.0) }
###		or
###		random = Random.new
###		random.poisson(10)
###################################################################

###################################################################
####		組み合わせ数と階乗の計算
###################################################################
def calc_casen(n, r) 
	# 組み合わせの数を計算　n個からr個を選ぶ組合せを計算
    return 0	if n < r
	r = n - r	if n - r < r
    return 1	if r == 0
    return n	if r == 1

	a = []
    1.step(r - 1){|i| a[i] = i + 2 } 
    3.step(n - r + 1){|i| 
        a[0] = i 
        1.step(r - 1){|j| a[j] += a[j-1] } 
    } 
    return a[r - 1] 
end

def eval_binomial_pmf(k, n, p)
	# Evaluates the binomial pmf.
	# Returns the probabily of k successes in n trials with probability p.
    return calc_casen(n, k) * p**k * (1.0 - p)**(n - k) 
end

# 階乗計算(reduce)
def factr(n)
	return 1  if n == 0
	return (1..n).reduce(&:*)
end

def binomial_coef(n, k)
    # n から k を取り出す組合せの計算
    return calc_casen(n, k).to_f
end

def log_binomial_coef(n, k)
    # Computes the log of the binomial coefficient.
	#  http://math.stackexchange.com/questions/64716/ 		
	# 			approximating-the-logarithm-of-the-binomial-coefficient
    return n * Math.log(n) - k * Math.log(k) - (n - k) * Math.log(n - k)
end
#######################################################################
#####	 		以下 Plot ユーティリティ
#######################################################################
def make_plot_2d(t, title: "title", x_label:"x_data", y_label: "y_label",
				 file_name: nil, key: "right", grid: nil, size: nil, xrange: nil,
				 yrange: nil, zrange: nil, xtics: nil, ytics: nil, ztics: nil,
				 logscale: nil)
	tt = "plot t[0]"
	t[1..-1].each_index{|i|  tt += ",t[#{i+1}]" }
    Numo.gnuplot do
   	   	set terminal: "gif"   if file_name
   	   	set output: file_name   if file_name
   	   	set key: key 	if key			# "right" or "left" or "outside" etc
   	   	set nokey: ""	if key.nil?
   	   	set grid: ""	if grid
   	   	set size: size	if size			# "square" or "ratio 0.5"
   	   	set "autoscale"
   	   	set xrange: xrange	if xrange	# "[min:max]"
   	   	set yrange: yrange	if yrange
   	   	set xtics: xtics	if xtics	# "start,incr,end" or "(0,1,2,4,8)" + rotate
   	   	set ytics: ytics	if ytics
   	   	if logscale
  	   		unset "logscale"
   	   		set "logscale x"	if logscale =~ /[xX]/
   	   		set "logscale y"	if logscale =~ /[yY]/
  		else
  	   		unset "logscale"
  	   	end
  	   	set title: title
		set xlabel: x_label
		set ylabel: y_label
		eval(tt)
	    sleep (5)
    end
end

def make_plot_3d(t, title: "title", x_label:"x_data", y_label: "y_label", z_label: "z_label",
				 file_name: nil, key: "right", grid: nil, size: nil, xrange: nil,
				 yrange: nil, zrange: nil, xtics: nil, ytics: nil, ztics: nil,
				 logscale: nil, contour: nil, map: nil)
	tt = "splot t[0]"
	t[1..-1].each_index{|i|  tt += ",t[#{i+1}]" }
    Numo.gnuplot do
   	   	set terminal: "gif"   if file_name
   	   	set output: file_name   if file_name
   	   	set key: key 	if key			# "right" or "left" or "outside" etc
   	   	set nokey: ""	if key.nil?
   	   	set grid: ""	if grid
   	   	set size: size	if size			# "square" or "ratio 0.5"
   	   	set "autoscale"
   	   	set xrange: xrange	if xrange	# "[min:max]"
   	   	set yrange: yrange	if yrange
   	   	set zrange: zrange	if zrange
   	   	set xtics: xtics	if xtics	# "start,incr,end" or "(0,1,2,4,8)" + rotate
   	   	set ytics: ytics	if ytics
   	   	set ztics: ztics	if ztics
   	   	if logscale
  	   		unset "logscale"
   	   		set "logscale x"	if logscale =~ /[xX]/
   	   		set "logscale y"	if logscale =~ /[yY]/
   	   		set "logscale z"	if logscale =~ /[zZ]/
  		else
  	   		unset "logscale"
  	   	end
  	   	set title: title
		set xlabel: x_label
		set ylabel: y_label
		set zlabel: z_label
		set "hidden3d"
		set "dgrid3d 100,100,4"
		if map then set "pm3d map" else unset "view" end
		if contour then set contour: contour  else unset "contour" end	# "base" or "surface" or "both"
		eval(tt)
	    sleep (5)
    end
end
