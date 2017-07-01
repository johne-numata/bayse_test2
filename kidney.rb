#! ruby -Ks
require './thinkbayes'
require "numo/gnuplot"
#require 'randomext'

INTERVAL = 245/365.0
MINSIZE = 0.2
MAXSIZE = 20
BUCKET_FACTOR = 10
$rdt_flag = 0
$xx

def simple_model()
    # Runs calculations based on a simple model.
    # time between discharge and diagnosis, in days
    interval = 3291.0

    # doubling time in linear measure is doubling time in volume * 3
    dt = 811.0 * 3

    # number of doublings since discharge
    doublings = interval / dt

    # how big was the tumor at time of discharge (diameter in cm)
    d1 = 15.5
    d0 = d1 / 2.0 ** doublings

    print 'interval (days) = ', interval, "\n"
    print 'interval (years) = ', interval / 365, "\n"
    print 'dt = ', dt, "\n"
    print 'doublings = ', doublings, "\n"
    print 'd1 = ', d1, "\n"
    print 'd0 = ', d0, "\n"

    # assume an initial linear measure of 0.1 cm
    d0 = 0.1
    d1 = 15.5

    # how many doublings would it take to get from d0 to d1
    doublings = Math.log(d1 / d0, 2)

    # what linear doubling time does that imply?
    dt = interval / doublings

    print 'doublings = ', doublings, "\n"
    print 'dt = ', dt, "\n"

    # compute the volumetric doubling time and RDT
    vdt = dt / 3
    rdt = 365 / vdt

    print 'vdt = ', vdt, "\n"
    print 'rdt = ', rdt, "\n"

    cdf = make_cdf()
    p = cdf.prob(rdt)
    print 'Prob{RDT > 2.4} = ', 1 - p, "\n\n"
end

def make_cdf()
	# Uses the data from Zhang et al. to construct a CDF.
    n = 53.0
    freqs = [0, 2, 31, 42, 48, 51, 52, 53]
    ps = freqs.collect{|freq| freq / n }
    xs = -1.5.step(6, 1.0).to_a

    cdf = Cdf.new(xs, ps)
    return cdf
end

def plot_cdf(cdf)
    # Plots the actual and fitted distributions.
	#   cdf: CDF object
    xs, ps = cdf.xs, cdf.ps
    cps = ps.collect{|p| 1 - p }

    # CCDF on logy scale: shows exponential behavior
	make_plot_2d([xs, cps, w:'lines', t:'bo-'], x_label: 'RDT',
					 y_label: 'CCDF (log scale)', title: 'kidney1', logscale: "y", grid: "on")

    # CDF, model and data
    mxs, mys = model_cdf()
	make_plot_2d([xs, cps, w:'lines', t:'model'], x_label: 'RDT',
					 y_label: 'CCDF (log scale)', logscale: "y", grid: "on")

	make_plot_2d([xs, ps, w:'lines', t:'data'], x_label: 'RDT (volume doublings per year)',
					 y_label: "CDF", title: 'kidney2', xtics: '-2, 1, 7', grid: "on")
end

def qq_plot(cdf, fit)
	# Makes a QQPlot of the values from actual and fitted distributions.
	#   cdf: actual Cdf of RDT
    #   fit: model
    xs = [-1.5, 5.5]

    xs, ps = cdf.xs, cdf.ps
    fs = ps.collect{|p| fit.value(p) }
	make_plot_2d([[xs, xs, w:'line', t:"45Åã"],[xs, fs, w:'lines', t:'Model']],
					 x_label: 'RDT (volume doublings per year)', y_label: 'Actual', title: 'kidney3', grid: "on")
end

def fit_cdf(cdf)
    # Fits a line to the log CCDF and returns the slope.
	#   cdf: Cdf of RDT
    xs, ps = cdf.xs, cdf.ps
    cps = ps.collect{|p| 1 - p }
    xs = xs[1..-2]
    lcps = cps[1..-2].collect{|p| Math.log(p) }
    a = xs.reg_line(lcps)
    return -a[:slope]
end

def correlated_generator(cdf, rho)
    # Generates a sequence of values from cdf with correlation.
	# Generates a correlated standard Gaussian series, then transforms to
    # values from cdf
	#    cdf: distribution to choose from
    #    rho: target coefficient of correlation
    def transform(x, cdf)
        # Maps from a Gaussian variate to a variate with the given CDF.
        p = gaussian_cdf(x)
        return cdf.value(p)
	end
    # for the first value, choose from a Gaussian and transform it
    r = Random.new()
	if $rdt_flag == 0
    	$rdt_flag = 1
    	return transform($xx = r.normal(0, 1), cdf)
    else 
    # for subsequent values, choose from the conditional distribution
    # based on the previous value
    	sigma = Math.sqrt(1 - rho**2)
    	return transform($xx = r.normal($xx * rho, sigma), cdf)
	end
end

def uncorrelated_generator(cdf, _rho = nil)
    # Generates a sequence of values from cdf with no correlation.
	# Ignores rho, which is accepted as a parameter to provide the
    # same interface as CorrelatedGenerator
	#   cdf: distribution to choose from
    #   rho: ignored
        return cdf.random()
end

def rdt_generator(cdf, rho)
    # Returns an iterator with n values from cdf and the given correlation.
	#    cdf: Cdf object
    #    rho: coefficient of correlation
    if rho == 0.0
        return uncorrelated_generator(cdf)
    else
        return correlated_generator(cdf, rho)
	end
end

def generate_rdt(pc, lam1, lam2)
    # Generate an RDT from a mixture of exponential distributions.
	# With prob pc, generate a negative value with param lam2;
    # otherwise generate a positive value with param lam1.
    return  rand < pc ? -(Random::DEFAULT.exponential(lam2)) : Random::DEFAULT.exponential(lam1)
end

def generate_sample(n, pc, lam1, lam2)
    # Generates a sample of RDTs.
	#    n: sample size
    #    pc: probablity of negative growth
    #    lam1: exponential parameter of positive growth
    #    lam2: exponential parameter of negative growth
	#    Returns: list of random variates
    xs = n.times.collect{generate_rdt(pc, lam1, lam2)}
    return xs
end

def generate_cdf(lam1 = 0.79, lam2 = 5.0, n = 1000, pc = 0.35)
    # Generates a sample of RDTs and returns its CDF.
	#    n: sample size
    #    pc: probablity of negative growth
    #    lam1: exponential parameter of positive growth
    #    lam2: exponential parameter of negative growth
	#    Returns: Cdf of generated sample
    xs = generate_sample(n, pc, lam1, lam2)
    cdf = make_cdf_from_list(xs)
    return cdf
end

def model_cdf(pc = 0.35, lam1 = 0.79, lam2 = 5.0)
    #    pc: probablity of negative growth
    #    lam1: exponential parameter of positive growth
    #    lam2: exponential parameter of negative growth
	#    Returns: list of xs, list of ys
    x1 = -2.step(0, 0.1).to_a
    y1 = x1.collect{|x| pc * (1 - eval_exponential_cdf(-x, lam2)) }
    x2 = 0.step(7, 0.1).to_a
    y2 = x2.collect{|x| pc + (1-pc) * eval_exponential_cdf(x, lam1) }
    return x1 + x2, y1 + y2
end

def bucket_to_cm(y, factor = BUCKET_FACTOR)
    # Computes the linear dimension for a given bucket.
	#   y: bucket number
    #   factor: multiplicitive factor from one bucket to the next
	#   Returns: linear dimension in cm
    return Math.exp(y.to_f / factor)
end

def cm_to_bucket(x, factor = BUCKET_FACTOR)
    # Computes the bucket for a given linear dimension.
	#   x: linear dimension in cm
    #   factor: multiplicitive factor from one bucket to the next
	#   Returns: float bucket number
	return -200 if x == 0
    return (factor * Math.log(x)).round
end

def diameter(volume, factor = 3 / Math::PI / 4, exp = 1 / 3.0)
    # Converts a volume to a diameter.
	#   d = 2r = 2 * (3/4/pi V)^1/3
    return 2 * (factor * volume) ** exp
end

def volume(diameter, factor = 4 * Math::PI / 3)
    # Converts a diameter to a volume.
	#   V = 4/3 pi (d/2)^3
    return factor * (diameter / 2.0)**3
end

class Cache
    # Records each observation point for each tumor.
	def initialize
		#   joint: map from (age, bucket) to frequency
        #   sequences: map from bucket to a list of sequences
        #   initial_rdt: sequence of (V0, rdt) pairs
        @joint = Joint.new()
        @sequences = {}
        @initial_rdt = []
	end

    def get_buckets()
        # Returns an iterator for the keys in the cache.
        return @sequences.keys
	end

    def get_sequence(bucket)
        # Looks up a bucket in the cache.
        return @sequences[bucket]
    end

    def conditional_cdf(bucket)
        # Forms the cdf of ages for a given bucket.
		#     bucket: int bucket number
        #     name: string
        pmf = @joint.conditional(0, 1, bucket)
        cdf = pmf.make_cdf()
        return cdf
    end

    def prob_older(cm, age)
        # Computes the probability of exceeding age, given size.
		#    cm: size in cm
        #    age: age in years
        bucket = cm_to_bucket(cm)
        cdf = conditional_cdf(bucket)
        p = cdf.prob(age)
        return 1 - p
    end

    def get_dist_age_size(size_thresh = MAXSIZE)
        # Gets the joint distribution of age and size.
		# Map from (age, log size in cm) to log freq
        #    Returns: new Pmf object
        joint = Joint.new()

        @joint.each{|val, freq|
            age, bucket = val
            cm = bucket_to_cm(bucket)
            next if cm > size_thresh
            log_cm = Math.log(cm, 10)
            joint.set([age, log_cm], Math.log(freq) * 10)
		}
        return joint
	end

    def add(age, seq, rdt)
        # Adds this observation point to the cache.
		#   age: age of the tumor in years
        #   seq: sequence of volumes
        #   rdt: RDT during this interval
        final = seq[-1]
        cm = diameter(final)
        bucket = cm_to_bucket(cm)
        @joint.incr([age, bucket])

        @sequences.include?(bucket) ? @sequences[bucket] << seq : @sequences[bucket] = [seq]

        initial = seq[-2]
        @initial_rdt << [initial, rdt]
	end

    def print
        # Prints the size (cm) for each bucket, and the number of sequences.
        get_buckets.sort.each{|bucket|
            ss = get_sequence(bucket)
            diameter = bucket_to_cm(bucket)
            print diameter, ss.size
        }
    end	

    def correlation
        # Computes the correlation between log volumes and rdts.
        vs, rdts = @initial_rdt.transpose
        lvs = vs.collect{|v| Math.log(v) }
#        return correlation.Corr(lvs, rdts) # Å©ëää÷åWêîÇÃåvéZ
        return lvs.corr(rdts)				# Å©ëää÷åWêîÇÃåvéZ
    end
end

class Calculator
    # Encapsulates the state of the computation.
    attr_reader  :cache
    def initialize
        @cache = Cache.new()
	end

    def make_sequences(n, rho, cdf)
        # Returns a list of sequences of volumes.
		#   n: number of sequences to make
        #   rho: serial correlation
        #   cdf: Cdf of rdts
		#   Returns: list of n sequences of volumes
        sequences = []
        $rdt_flag = 0			# ÉäÉZÉbÉgå„íˆïœçX
        n.times{|i|
            rdt_seq = Object.method(:rdt_generator)
            seq = make_sequence(rdt_seq, cdf, rho)
            sequences << seq

            print "i= ", i, "\n"  if i % 100 == 0
		}
        return sequences
    end

    def make_sequence(rdt_seq, cdf, rho, v0 = 0.01, interval = INTERVAL, vmax = volume(MAXSIZE))
        # Simulate the growth of a tumor.
		#   rdt_seq: sequence of rdts
        #   v0: initial volume in mL (cm^3)
        #   interval: timestep in years
        #   vmax: volume to stop at
		#   Returns: sequence of volumes
        seq = [v0]
        age = 0
        while true do
        	rdt = rdt_seq.call(cdf, rho)
            age += interval
            final, seq = extend_sequence(age, seq, rdt, interval)
            break if final > vmax
        end
        return seq
	end

    def extend_sequence(age, seq, rdt, interval)
        # Generates a new random value and adds it to the end of seq.
		#    Side-effect: adds sub-sequences to the cache.
		#    age: age of tumor at the end of this interval
        #    seq: sequence of values so far
        #    rdt: reciprocal doubling time in doublings per year
        #    interval: timestep in years
		#    Returns: final volume, extended sequence
        initial = seq[-1]
        doublings = rdt * interval
        final = initial * 2**doublings
        new_seq = seq + [final]
        @cache.add(age, new_seq, rdt)

        return final, new_seq
	end

    def plot_buckets()
        # Plots the set of sequences that ended in a given bucket."""
        # 2.01, 4.95 cm, 9.97 cm
        buckets = [7.0, 16.0, 23.0]
        buckets.each{|a|
        	sequences = @cache.get_sequence(a.to_i)
        	plot_sequences(sequences)
        }
	end

    def plot_joint_dist()
        # Makes a pcolor plot of the age-size joint distribution.
        joint = @cache.get_dist_age_size()
		dd = joint.collect{|v| v.flatten}.transpose
		Numo.gnuplot do
			set title:'kidney8'
			set xlabel:'ages'
  			set ylabel:'diameter (cm, log scale)'
			splot dd[0], dd[1], dd[2]
		sleep (10)
		end
#        thinkplot.Contour(joint, contour=False, pcolor=True)

#        thinkplot.Save(root='kidney8',
#                    formats=FORMATS,
#                    axis=[0, 41, -0.7, 1.31],
#                    yticks=MakeLogTicks([0.2, 0.5, 1, 2, 5, 10, 20]),
#                    xlabel='ages',
#                    ylabel='diameter (cm, log scale)')
	end

    def plot_conditional_cdfs()
        # Plots the cdf of ages for each bucket.
        buckets = [7.0, 16.0, 23.0, 27.0]
        # 2.01, 4.95 cm, 9.97 cm, 14.879 cm
        names = ['2 cm', '5 cm', '10 cm', '15 cm']
        t = []
        buckets.zip(names).each{|bucket, name|
            cdf = @cache.conditional_cdf(bucket); cdf.name = name
            t << cdf.render
		}
		make_plot_2d(t, x_label: 'tumor age (years)', y_label: 'CDF',
						 title: 'Distribution of age for several diameters', grid: "on")
	end

    def plot_credible_intervals(xscale = 'linear')
        # Plots the confidence interval for each bucket.
        xs = []
        ts = []
        percentiles = [95, 75, 50, 25, 5]
        min_size = 0.3
		labels = ['95th', '75th', '50th', '25th', '5th']

        # loop through the buckets, accumulate
        # xs: sequence of sizes in cm
        # ts: sequence of percentile tuples
        @cache.get_buckets.sort!.each{|bucket|
            cm = bucket_to_cm(bucket)
            next  if cm < min_size or cm > 20.0
            xs << cm
            cdf = @cache.conditional_cdf(bucket)      
            ps = percentiles.collect{|p| cdf.percentile(p)}
            ts << ps
		}

        # transpose the ts so we have sequences for each percentile rank
        yys = ts.transpose
		data = []
		yys.each_index{|i|
            # plot the data points
            data << [xs, yys[i], w: 'line', t: labels[i]]

            # plot the fit lines
            fxs = [min_size, 20.0]
            fys = fit_line(xs, yys[i], fxs)

            data << [fxs, fys, w:'line']
		}
		make_plot_2d(data, x_label: 'tumor age (years)', y_label: 'diameter (cm, log scale)',
						 title: 'Credible interval for age vs diameter', grid: "on")
#                       xticks=MakeTicks([0.5, 1, 2, 5, 10, 20]),
#                       axis=[0.25, 35, 0, 45],
	end
end

def plot_sequences(sequences)
    # Plots linear measurement vs time.
	#   sequences: list of sequences of volumes
#    options = dict(color='gray', linewidth=1, linestyle='dashed')
#    thinkplot.Plot([0, 40], [10, 10], **options)
	t = []
    sequences.each{|seq|
		n = seq.size
        age = n * INTERVAL
        ts = 0.step(age, age / (n - 1)).to_a
		xs = seq.collect{|v| diameter(v)}
      	t << [ts, xs, {w:'lines'}]
	}
	tt = "plot t[0]"
	t[1..-1].each_index{|i|  tt += ",t[#{i+1}]" }
    Numo.gnuplot do
   	   	set title: 'kidney4'
		set xlabel: 'tumor age (years)'
		set ylabel: 'diameter (cm, log scale)'
		set 'logscale y'
       	set xtics: '0, 5, 40'
#		MakeTicks([0.2, 0.5, 1, 2, 5, 10, 20])
#       axis=[0, 40, MINSIZE, 20],
		eval(tt)
	    sleep (5)
    end
end

def print_ci(fp, cm, ps)
    # Writes a line in the LaTeX table.
	#    fp: file pointer
    #    cm: diameter in cm
    #    ts: tuples of percentiles
    fp.write('%0.1f' % round(cm, 1))
    for p in reversed(ps)
        fp.write(' & %0.1f ' % round(p, 1))
    end
    fp.write(r'\\' '\n')
end

def print_table(fp, xs, ts)
    # Writes the data in a LaTeX table.
	#    fp: file pointer
    #    xs: diameters in cm
    #    ts: sequence of tuples of percentiles
    fp.write(r'\begin{tabular}{|r||r|r|r|r|r|}' '\n')
    fp.write(r'\hline' '\n')
    fp.write(r'Diameter   & \multicolumn{5}{c|}{Percentiles of age} \\' '\n')
    fp.write(r'(cm)   & 5th & 25th & 50th & 75th & 95th \\' '\n')
    fp.write(r'\hline' '\n')

    for i, (cm, ps) in enumerate(zip(xs, ts))
        #print cm, ps
        if i % 3 == 0
            print_ci(fp, cm, ps)
        end
	end
    fp.write(r'\hline' '\n')
    fp.write(r'\end{tabular}' '\n')
end

def fit_line(xs, ys, fxs)
    # Fits a line to the xs and ys, and returns fitted values for fxs.
	# Applies a log transform to the xs.
	#    xs: diameter in cm
    #    ys: age in years
    #    fxs: diameter in cm
    lxs = xs.collect{|x| Math.log(x)}
    ans = lxs.reg_line(ys)
    # res = correlation.Residuals(lxs, ys, inter, slope)
    # r2 = correlation.CoefDetermination(ys, res)

    lfxs = fxs.collect{|x| Math.log(x)}
    fys = lfxs.collect{|x| ans[:intercept] + ans[:slope] * x}
    return fys
end

def test_correlation(cdf)
    # Tests the correlated generator.
	# Makes sure that the sequence has the right distribution and correlation.
    n = 10000
    rho = 0.4

     xs = n.times.collect{ correlated_generator(cdf, rho)}
#    rho2 = correlation.SerialCorr(xs)  # ånóÒëää÷åWêîÇÃåvéZ Pythonî≈
    rho2 = xs.serial_corr				# ånóÒëää÷åWêîÇÃåvéZÅAarrayÇÃägí£ 
    print 'rho = ', rho, "  ", 'rho2 = ', rho2, "\n"
    cdf2 = make_cdf_from_list(xs)

#    make_plot_cdfs([["ORG", cdf], ["NEW", cdf2]])
end

################################################
###    main ÉãÅ[É`Éì
################################################
for size in [1, 5, 10]
    bucket = cm_to_bucket(size)
    print 'Size, bucket = ', size, "  ", bucket, "\n"
end
simple_model()

cdf = make_cdf()
lam1 = fit_cdf(cdf)
fit = generate_cdf(lam1)
test_correlation(fit)

plot_cdf(cdf)
qq_plot(cdf, fit)
calc = Calculator.new()
rho = 0.0
sequences = calc.make_sequences(100, rho, fit)
#plot_sequences(sequences)
#rho = 0.6
#sequences = calc.make_sequences(100, rho, fit)
#plot_sequences(sequences)

calc.plot_buckets()
rho = 0.0
calc.make_sequences(1900, rho, fit)
print 'V0-RDT correlation = ', calc.cache.correlation(), "\n"
print '15.5 Probability age > 8 year = ', calc.cache.prob_older(15.5, 8), "\n"
print '6.0 Probability age > 8 year = ', calc.cache.prob_older(6.0, 8), "\n"

calc.plot_conditional_cdfs()
calc.plot_credible_intervals(xscale = 'log')
calc.plot_joint_dist()
