#######################
##    Gaussian KDE Class
##    æŠ¸‚¦‚¸1ŸŒ³ŒÀ’è
##    2015-04-19
#######################
require 'randomext'
class GaussianKde
#    Kernel density estimation is a way to estimate the probability density
#    function (PDF) of a random variable in a non-parametric way.
#    `gaussian_kde` works for both uni-variate and multi-variate data.   It
#    includes automatic bandwidth determination.  The estimation works best for
#    a unimodal distribution; bimodal or multi-modal distributions tend to be
#    oversmoothed.
#
#    Parameters
#    ----------
#    dataset : array
#    bw_method : str, scalar or callable, optional
#        The method used to calculate the estimator bandwidth.  This can be
#        'scott', 'silverman', a scalar constant.  If a scalar, this will be
#		 used directly as `@band_width`.  
#        If None (default), 'silverman' is used.  See Notes for more details.
#
#    Attributes
#    ----------
#    @dataset : array
#        The dataset with which `gaussian_kde` was initialized.
#    @band_width : float
#        The bandwidth factor.
#
#    Methods
#    -------
#    kde.evaluate(points) : ndarray
#        Evaluate the estimated pdf on a provided set of points.
#    kde.integrate_gaussian(mean, cov) : float
#        Multiply pdf with a specified Gaussian and integrate over the whole
#        domain.
#    kde.integrate_box(low, high) : float
#        Integrate pdf (1D only) between two bounds.
#    kde.resample(size=None) : array
#        Randomly sample a dataset from the estimated pdf.
#    kde.set_bandwidth(bw_method='scott') : None
#        Computes the bandwidth, i.e. the coefficient that multiplies the data
#        covariance matrix to obtain the kernel covariance matrix.
#        .. versionadded:: 0.11.0
#

    def initialize(dataset, bw_method="silverman")
		begin 					#ˆêŸŒ³‚Ìê‡‚Í“ñŸŒ³‚Ö•ÏŠ·
			dataset.dig(1,1)
		rescue
#			puts "error \n"
			dataset = dataset.zip([].fill(0..dataset.size - 1){"1"})
		end
        @dataset = dataset
        @size = dataset.size
        if @size <= 1
            raise "`dataset` input should have multiple elements."
		end
        set_bandwidth(@bw_method = bw_method)
	end

    def evaluate(points)
        number = points.size
        result = [].fill(0.0, 0..number -1)

        for i in 0..number - 1
            diff = @dataset.map{|item| (item[0].to_f - points[i].to_f) / @band_width}
            result[i] += diff.inject(0.0){|sum, item| sum + Math.exp(item ** 2 / -2.0)}
		end
		sum = result.sum
        result.map!{|item| item / sum}

        return result
	end

#    def integrate_gaussian(self, mean, cov):
#        """
#        Multiply estimated density by a multivariate Gaussian and integrate
#        over the whole space.
#
#        Parameters
#        ----------
#        mean : aray_like
#            A 1-D array, specifying the mean of the Gaussian.
#        cov : array_like
#            A 2-D array, specifying the covariance matrix of the Gaussian.
#
#        Returns
 #       -------
#        result : scalar
#            The value of the integral.
#        Raises
#        ------
#        ValueError :
#            If the mean or covariance of the input Gaussian differs from
#            the KDE's dimensionality.
#
#        """
#        mean = atleast_1d(squeeze(mean))
#        cov = atleast_2d(cov)
#
#        if mean.shape != (self.d,):
#            raise ValueError("mean does not have dimension %s" % self.d)
 #       if cov.shape != (self.d, self.d):
#            raise ValueError("covariance does not have dimension %s" % self.d)
#
#        # make mean a column vector
#        mean = mean[:, newaxis]
#
#        sum_cov = self.covariance + cov
#
#        diff = self.dataset - mean
 #       tdiff = dot(linalg.inv(sum_cov), diff)
#
#        energies = sum(diff * tdiff, axis=0) / 2.0
#        result = sum(exp(-energies), axis=0) / sqrt(linalg.det(2 * pi *
#                                                        sum_cov)) / self.n
#
#        return result
#
#    def integrate_box(low, high)
#		###  Computes the integral of a 1D pdf between two bounds.
#        stdev = ravel(sqrt(self.covariance))[0]
#
#        normalized_low = ravel((low - self.dataset) / stdev)
#        normalized_high = ravel((high - self.dataset) / stdev)
#
#        value = np.mean(special.ndtr(normalized_high) -
#                        special.ndtr(normalized_low))
#        return value
#		end

    def resample(request = nil)
        ###  Randomly sample a dataset from the estimated pdf.
        if request == nil
            request = @size
		end
		result = []
        mu = 0
        sigma = @band_width
        # ³‹K•ª•z‚É]‚¤—”‚ğ‘«‚·
		request.times{|i| result[i] = @dataset[rand(@size)][0].to_f + Random::DEFAULT.normal(mu,sigma) }
		
        return result
	end

    def scotts_factor()
    	cdf = make_cdf_from_items(@dataset)
    	iqr = (cdf.percentile(75) - cdf.percentile(25)) / 1.34
		sum = 0; var = 0
		@dataset.each{|a, b| sum += a.to_f * b.to_f}
		mean = sum / @dataset.size
		var = @dataset.reduce(0.0){|sum, b| sum + (b[0].to_f - mean) ** 2 * b[1].to_f}
		sd = Math.sqrt(var / (@dataset.size - 1))
        return [iqr, sd].min * 1.06 * @size**(-1.0/5)
	end

    def silverman_factor()
    	cdf = make_cdf_from_items(@dataset)
    	iqr = (cdf.percentile(75) - cdf.percentile(25)) / 1.34
		sum = 0; var = 0											# ˆÈ~sd‚ÌŒvZ
		@dataset.each{|a, b| sum += a.to_f * b.to_f}
		mean = sum / @dataset.size
		var = @dataset.reduce(0.0){|sum, b| sum + (b[0].to_f - mean) ** 2 * b[1].to_f}
		sd = Math.sqrt(var / (@dataset.size - 1))
        return [iqr, sd].min * 0.9 * @size**(-1.0/5)
	end

    def set_bandwidth(bw_method = 'silverman')
		###  Compute the estimator bandwidth with given method.
        if bw_method == 'scott'
            @band_width = scotts_factor()
        elsif bw_method == 'silverman'
            @band_width = silverman_factor()
        elsif !(bw_method.nil?)
            @band_width = bw_method.to_f
        else
            raise "`bw_method` should be 'scott', 'silverman', a scalar "
        end
	end
end	
