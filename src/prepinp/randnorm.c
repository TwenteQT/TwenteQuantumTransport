#include <errno.h>
#include <math.h>

/* Coefficients in rational approximations. */
static const double randnc_a[] =
{
	-3.969683028665376e+01,
	 2.209460984245205e+02,
	-2.759285104469687e+02,
	 1.383577518672690e+02,
	-3.066479806614716e+01,
	 2.506628277459239e+00
};

static const double randnc_b[] =
{
	-5.447609879822406e+01,
	 1.615858368580409e+02,
	-1.556989798598866e+02,
	 6.680131188771972e+01,
	-1.328068155288572e+01
};

static const double randnc_c[] =
{
	-7.784894002430293e-03,
	-3.223964580411365e-01,
	-2.400758277161838e+00,
	-2.549732539343734e+00,
	 4.374664141464968e+00,
	 2.938163982698783e+00
};

static const double randnc_d[] =
{
	7.784695709041462e-03,
	3.224671290700398e-01,
	2.445134137142996e+00,
	3.754408661907416e+00
};

#define randnc_LOW 0.02425
#define randnc_HIGH 0.97575

double ltqnorm(double p)
{
	double q, r;

	errno = 0;

	if (p < 0 || p > 1)
	{
		errno = EDOM;
		return 0.0;
	}
	else if (p == 0)
	{
		errno = ERANGE;
		return -HUGE_VAL /* minus "infinity" */;
	}
	else if (p == 1)
	{
		errno = ERANGE;
		return HUGE_VAL /* "infinity" */;
	}
	else if (p < randnc_LOW)
	{
		/* Rational approximation for lower region */
		q = sqrt(-2*log(p));
		return (((((randnc_c[0]*q+randnc_c[1])*q+randnc_c[2])*q+randnc_c[3])*q+randnc_c[4])*q+randnc_c[5]) /
			((((randnc_d[0]*q+randnc_d[1])*q+randnc_d[2])*q+randnc_d[3])*q+1);
	}
	else if (p > randnc_HIGH)
	{
		/* Rational approximation for upper region */
		q  = sqrt(-2*log(1-p));
		return -(((((randnc_c[0]*q+randnc_c[1])*q+randnc_c[2])*q+randnc_c[3])*q+randnc_c[4])*q+randnc_c[5]) /
			((((randnc_d[0]*q+randnc_d[1])*q+randnc_d[2])*q+randnc_d[3])*q+1);
	}
	else
	{
		/* Rational approximation for central region */
    		q = p - 0.5;
    		r = q*q;
		return (((((randnc_a[0]*r+randnc_a[1])*r+randnc_a[2])*r+randnc_a[3])*r+randnc_a[4])*r+randnc_a[5])*q /
			(((((randnc_b[0]*r+randnc_b[1])*r+randnc_b[2])*r+randnc_b[3])*r+randnc_b[4])*r+1);
	}
}
