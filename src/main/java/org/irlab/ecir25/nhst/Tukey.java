package org.irlab.ecir25.nhst;

import org.apache.commons.math3.distribution.NormalDistribution;

import static java.lang.Math.*;

// Implementation of the Tukey CDF to be used in the TukeyHSD test.
// All this code is copied from R sources.
public class Tukey {

  public static final double M_PI = 3.141592653589793238462643383279502884197169399375;
  public static final double M_1_SQRT_2PI = 0.398942280401432677939946059934;
  public static final double M_LN2 = 0.693147180559945309417232121458;
  public static final double M_LN_SQRT_2PI = 0.91893853320467274178032973640562;
  public static final double M_LN_SQRT_PId2 = 0.225791352644727432363097614947441;

  private static double wprob(double w, double rr, double cc) {
    final int nleg = 12, ihalf = 6;

		/* looks like this is suboptimal for double precision.
	       (see how C1-C3 are used) <MM>
		 */
    /* const double iMax  = 1.; not used if = 1*/
    final double C1 = -30.;
    final double C2 = -50.;
    final double C3 = 60.;
    final double bb = 8.;
    final double wlar = 3.;
    final double wincr1 = 2.;
    final double wincr2 = 3.;
    final double[] xleg = { 0.981560634246719250690549090149,
                            0.904117256370474856678465866119,
                            0.769902674194304687036893833213,
                            0.587317954286617447296702418941,
                            0.367831498998180193752691536644,
                            0.125233408511468915472441369464 };
    final double[] aleg = { 0.047175336386511827194615961485,
                            0.106939325995318430960254718194,
                            0.160078328543346226334652529543,
                            0.203167426723065921749064455810,
                            0.233492536538354808760849898925,
                            0.249147045813402785000562436043 };
    double a, ac, pr_w, b, binc, c, cc1, pminus, pplus, qexpo, qsqz, rinsum, wi, wincr, xx;
    /* long */
    double blb, bub, einsum, elsum;
    int j, jj;


    qsqz = w * 0.5;

    /* if w >= 16 then the integral lower bound (occurs for c=20) */
    /* is 0.99999999999995 so return a value of 1. */

    if (qsqz >= bb) return 1.0;

    /* find (f(w/2) - 1) ^ cc */
    /* (first term in integral of hartley's form). */

    pr_w = 2 * new NormalDistribution(0, 1).cumulativeProbability(qsqz) - 1.; /* erf(qsqz / M_SQRT2) */
    /* if pr_w ^ cc < 2e-22 then set pr_w = 0 */
    if (pr_w >= exp(C2 / cc)) pr_w = pow(pr_w, cc);
    else pr_w = 0.0;

    /* if w is large then the second component of the */
    /* integral is small, so fewer intervals are needed. */

    if (w > wlar) wincr = wincr1;
    else wincr = wincr2;

    /* find the integral of second term of hartley's form */
    /* for the integral of the range for equal-length */
    /* intervals using legendre quadrature.  limits of */
    /* integration are from (w/2, 8).  two or three */
    /* equal-length intervals are used. */

    /* blb and bub are lower and upper limits of integration. */

    blb = qsqz;
    binc = (bb - qsqz) / wincr;
    bub = blb + binc;
    einsum = 0;

    /* integrate over each interval */

    cc1 = cc - 1.0;
    for (wi = 1; wi <= wincr; wi++) {
      elsum = 0;
      a = 0.5 * (bub + blb);


      /* legendre quadrature with order = nleg */

      b = 0.5 * (bub - blb);

      for (jj = 1; jj <= nleg; jj++) {
        if (ihalf < jj) {
          j = (nleg - jj) + 1;
          xx = xleg[j - 1];
        } else {
          j = jj;
          xx = -xleg[j - 1];
        }
        c = b * xx;
        ac = a + c;

        /* if exp(-qexpo/2) < 9e-14, */
        /* then doesn't contribute to integral */

        qexpo = ac * ac;
        if (qexpo > C3) break;


        pplus = 2 * new NormalDistribution(0, 1).cumulativeProbability(ac);
        pminus = 2 * new NormalDistribution(w, 1).cumulativeProbability(ac);

        /* if rinsum ^ (cc-1) < 9e-14, */
        /* then doesn't contribute to integral */

        rinsum = (pplus * 0.5) - (pminus * 0.5);
        if (rinsum >= exp(C1 / cc1)) {
          rinsum = (aleg[j - 1] * exp(-(0.5 * qexpo))) * pow(rinsum, cc1);

          elsum += rinsum;
        }
      }
      elsum *= (((2.0 * b) * cc) * M_1_SQRT_2PI);
      einsum += elsum;
      blb = bub;
      bub += binc;
    }

    /* if pr_w ^ rr < 9e-14, then return 0 */
    pr_w = einsum + pr_w;
    if (pr_w <= exp(C1 / rr)) return 0.;

    pr_w = pow(pr_w, rr);
    if (pr_w >= 1.)/* 1 was iMax was eps */ return 1.;
    return pr_w;
  } /* wprob() */

  public static double cumulative(double q, double rr, double cc, double df, boolean lower_tail, boolean log_p) {
    final int nlegq = 16, ihalfq = 8;

    /*  const double eps = 1.0; not used if = 1 */
    final double eps1 = -30.0;
    final double eps2 = 1.0e-14;
    final double dhaf = 100.0;
    final double dquar = 800.0;
    final double deigh = 5000.0;
    final double dlarg = 25000.0;
    final double ulen1 = 1.0;
    final double ulen2 = 0.5;
    final double ulen3 = 0.25;
    final double ulen4 = 0.125;
    final double[] xlegq = { 0.989400934991649932596154173450,
                             0.944575023073232576077988415535,
                             0.865631202387831743880467897712,
                             0.755404408355003033895101194847,
                             0.617876244402643748446671764049,
                             0.458016777657227386342419442984,
                             0.281603550779258913230460501460,
                             0.950125098376374401853193354250e-1 };
    final double[] alegq = { 0.271524594117540948517805724560e-1,
                             0.622535239386478928628438369944e-1,
                             0.951585116824927848099251076022e-1,
                             0.124628971255533872052476282192,
                             0.149595988816576732081501730547,
                             0.169156519395002538189312079030,
                             0.182603415044923588866763667969,
                             0.189450610455068496285396723208 };
    double ans, f2, f21, f2lf, ff4, otsum = 0, qsqz, rotsum, t1, twa1, ulen, wprb;
    int i, j, jj;

    if (Double.isInfinite(q) || Double.isInfinite(rr) || Double.isInfinite(cc) || Double.isInfinite(df))
      return Double.NaN;

    if (q <= 0) return (lower_tail ? (log_p ? Double.NEGATIVE_INFINITY : 0.) : (log_p ? 0. : 1.));

    /* df must be > 1 */
    /* there must be at least two values */

    if (df < 2 || rr < 1 || cc < 2) return Double.NaN;

    if (Double.isInfinite(q)) return (lower_tail ? (log_p ? 0. : 1.) : (log_p ? Double.NEGATIVE_INFINITY : 0.));

    if (df > dlarg) {
      //return R_DT_val(wprob(q, rr, cc));
      double x = wprob(q, rr, cc);
      return (lower_tail ? (log_p ? log(x) : (x)) : (log_p ? log1p(-(x)) : (0.5 - (x) + 0.5)));
    }

    /* calculate leading constant */

    f2 = df * 0.5;
    /* lgammafn(u) = log(gamma(u)) */
    f2lf = ((f2 * log(df)) - (df * M_LN2)) - lgammafn(f2);
    f21 = f2 - 1.0;

    /* integral is divided into unit, half-unit, quarter-unit, or */
    /* eighth-unit length intervals depending on the value of the */
    /* degrees of freedom. */

    ff4 = df * 0.25;
    if (df <= dhaf) ulen = ulen1;
    else if (df <= dquar) ulen = ulen2;
    else if (df <= deigh) ulen = ulen3;
    else ulen = ulen4;

    f2lf += log(ulen);

    /* integrate over each subinterval */

    ans = 0.0;

    for (i = 1; i <= 50; i++) {
      otsum = 0.0;

      /* legendre quadrature with order = nlegq */
      /* nodes (stored in xlegq) are symmetric around zero. */

      twa1 = (2 * i - 1) * ulen;

      for (jj = 1; jj <= nlegq; jj++) {
        if (ihalfq < jj) {
          j = jj - ihalfq - 1;
          t1 = (f2lf + (f21 * log(twa1 + (xlegq[j] * ulen)))) - (((xlegq[j] * ulen) + twa1) * ff4);
        } else {
          j = jj - 1;
          t1 = (f2lf + (f21 * log(twa1 - (xlegq[j] * ulen)))) + (((xlegq[j] * ulen) - twa1) * ff4);

        }

        /* if exp(t1) < 9e-14, then doesn't contribute to integral */
        if (t1 >= eps1) {
          if (ihalfq < jj) {
            qsqz = q * sqrt(((xlegq[j] * ulen) + twa1) * 0.5);
          } else {
            qsqz = q * sqrt(((-(xlegq[j] * ulen)) + twa1) * 0.5);
          }

          /* call wprob to find integral of range portion */

          wprb = wprob(qsqz, rr, cc);
          rotsum = (wprb * alegq[j]) * exp(t1);
          otsum += rotsum;
        }
        /* end legendre integral for interval i */
        /* L200: */
      }

      /* if integral for interval i < 1e-14, then stop.
       * However, in order to avoid small area under left tail,
       * at least  1 / ulen  intervals are calculated.
       */
      if (i * ulen >= 1.0 && otsum <= eps2) break;

      /* end of interval i */
      /* L330: */

      ans += otsum;
    }

    if (otsum > eps2) { /* not converged */
      //ML_ERROR(ME_PRECISION, "ptukey");
      System.err.println("Precision error at Tukey.cumulative");
    }
    if (ans > 1.) ans = 1.;
    //return R_DT_val(ans);
    return (lower_tail ? (log_p ? log(ans) : (ans)) : (log_p ? log1p(-(ans)) : (0.5 - (ans) + 0.5)));
  }

  private static final double lgammafn_sign(double x, int[] sgn) {
    final double xmax = 2.5327372760800758e+305;
    final double dxrel = 1.490116119384765696e-8;
    double ans, y, sinpiy;

    if (sgn != null) sgn[0] = 1;
    if (x < 0 && (trunc(-x) % 2.) == 0) if (sgn != null) sgn[0] = -1;

    if (x <= 0 && x == floor(x)) { /* Negative integer argument */
      //ML_ERROR(ME_RANGE, "lgamma");
      return Double.POSITIVE_INFINITY;/* +Inf, since lgamma(x) = log|gamma(x)| */
    }

    y = abs(x);

    if (y < 1e-306) return -log(x); // denormalized range, R change
    if (y <= 10) return log(abs(gammafn(x)));
    if (y > xmax) {
      //ML_ERROR(ME_RANGE, "lgamma");
      return Double.POSITIVE_INFINITY;
    }

    if (x > 0) { /* i.e. y = x > 10 */
      if (x > 1e17) return (x * (log(x) - 1.));
      if (x > 4934720.) return (M_LN_SQRT_2PI + (x - 0.5) * log(x) - x);
      return M_LN_SQRT_2PI + (x - 0.5) * log(x) - x + lgammacor(x);
    }
    /* else: x < -10; y = -x */
    sinpiy = abs(sin(M_PI * y));

    if (sinpiy == 0) { /* Negative integer argument ===
	    			  Now UNNECESSARY: caught above */
      //MATHLIB_WARNING(" ** should NEVER happen! *** [lgamma.c: Neg.int, y=%g]\n",y);
      return Double.NaN;
    }

    ans = M_LN_SQRT_PId2 + (x - 0.5) * log(y) - x - log(sinpiy) - lgammacor(y);

    if (abs((x - trunc(x - 0.5)) * ans / x) < dxrel) {
      /* The answer is less than half precision because
       * the argument is too near a negative integer. */
      //ML_ERROR(ME_PRECISION, "lgamma");
      System.err.println("lgamma precision error!");
    }
    return ans;
  }

  private static double trunc(double x) {
    return x >= 0 ? floor(x) : ceil(x);
  }


  private static double lgammacor(double x) {
    final double[] algmcs = { +.1666389480451863247205729650822e+0,
                              -.1384948176067563840732986059135e-4,
                              +.9810825646924729426157171547487e-8,
                              -.1809129475572494194263306266719e-10,
                              +.6221098041892605227126015543416e-13,
                              -.3399615005417721944303330599666e-15,
                              +.2683181998482698748957538846666e-17,
                              -.2868042435334643284144622399999e-19,
                              +.3962837061046434803679306666666e-21,
                              -.6831888753985766870111999999999e-23,
                              +.1429227355942498147573333333333e-24,
                              -.3547598158101070547199999999999e-26,
                              +.1025680058010470912000000000000e-27,
                              -.3401102254316748799999999999999e-29,
                              +.1276642195630062933333333333333e-30 };
    int nalgm = 5;
    double xbig = 94906265.62425156;
    double xmax = 3.745194030963158e306;
    double tmp;

    if (x < 10) return Double.NaN;
    if (x >= xmax) return 1 / (x * 12); // Underflow
    if (x < xbig) {
      tmp = 10 / x;
      return chebyshev_eval(tmp * tmp * 2 - 1, algmcs, nalgm) / x;
    }
    return 1 / (x * 12);
  }

  private static double chebyshev_eval(double x, double[] a, int n) {
    if (n < 1 || n > 1000 || x < -1.1 || x > 1.1) return Double.NaN;
    double twox = x * 2, b2 = 0, b1 = 0, b0 = 0;
    for (int i = 1; i <= n; i++) {
      b2 = b1;
      b1 = b0;
      b0 = twox * b1 - b2 + a[n - i];
    }
    return (b0 - b2) * 0.5;
  }

  private static double gammafn(double x) {
    final double[] gamcs = { +.8571195590989331421920062399942e-2,
                             +.4415381324841006757191315771652e-2,
                             +.5685043681599363378632664588789e-1,
                             -.4219835396418560501012500186624e-2,
                             +.1326808181212460220584006796352e-2,
                             -.1893024529798880432523947023886e-3,
                             +.3606925327441245256578082217225e-4,
                             -.6056761904460864218485548290365e-5,
                             +.1055829546302283344731823509093e-5,
                             -.1811967365542384048291855891166e-6,
                             +.3117724964715322277790254593169e-7,
                             -.5354219639019687140874081024347e-8,
                             +.9193275519859588946887786825940e-9,
                             -.1577941280288339761767423273953e-9,
                             +.2707980622934954543266540433089e-10,
                             -.4646818653825730144081661058933e-11,
                             +.7973350192007419656460767175359e-12,
                             -.1368078209830916025799499172309e-12,
                             +.2347319486563800657233471771688e-13,
                             -.4027432614949066932766570534699e-14,
                             +.6910051747372100912138336975257e-15,
                             -.1185584500221992907052387126192e-15,
                             +.2034148542496373955201026051932e-16,
                             -.3490054341717405849274012949108e-17,
                             +.5987993856485305567135051066026e-18,
                             -.1027378057872228074490069778431e-18,
                             +.1762702816060529824942759660748e-19,
                             -.3024320653735306260958772112042e-20,
                             +.5188914660218397839717833550506e-21,
                             -.8902770842456576692449251601066e-22,
                             +.1527474068493342602274596891306e-22,
                             -.2620731256187362900257328332799e-23,
                             +.4496464047830538670331046570666e-24,
                             -.7714712731336877911703901525333e-25,
                             +.1323635453126044036486572714666e-25,
                             -.2270999412942928816702313813333e-26,
                             +.3896418998003991449320816639999e-27,
                             -.6685198115125953327792127999999e-28,
                             +.1146998663140024384347613866666e-28,
                             -.1967938586345134677295103999999e-29,
                             +.3376448816585338090334890666666e-30,
                             -.5793070335782135784625493333333e-31 };
    int i, n;
    double y, value;
		/*
		int ngam = 0;
		double xmin = 0., xmax = 0., xsml = 0., dxrel = 0.;
		if (ngam == 0) {
			ngam = chebyshev_init(gamcs, 42, DBL_EPSILON/20);
			xmin=-170.5674972726612; xmax=171.61447887182298;
			xsml = exp(max(log(Double.MIN_VALUE), -log(Double.MAX_VALUE))+0.01);
			dxrel = sqrt(DBL_EPSILON);
		}
		/*/
    int ngam = 22;
    double xmin = -170.5674972726612, xmax = 171.61447887182298, xsml = 2.2474362225598545e-308, dxrel = 1.490116119384765696e-8;
    //*/

    if (Double.isNaN(x)) return x;
    // If the argument is exactly zero or a negative integer then return NaN.
    if (x == 0 || (x < 0 && x == (long) x)) return Double.NaN;
    y = abs(x);

    if (y <= 10) {
      /* Compute gamma(x) for -10 <= x <= 10. */
      /* Reduce the interval and find gamma(1 + y) for */
      /* 0 <= y < 1 first of all. */
      n = (int) x;
      if (x < 0) --n;
      y = x - n;/* n = floor(x)  ==>	y in [ 0, 1 ) */
      --n;
      value = chebyshev_eval(y * 2 - 1, gamcs, ngam) + .9375;
      if (n == 0) return value;/* x = 1.dddd = 1+y */

      if (n < 0) {
        /* compute gamma(x) for -10 <= x < 1 */

        /* The answer is less than half precision */
        /* because x too near a negative integer. */
        /*!* 	    if (x < -0.5 && fabs(x - (int)(x - 0.5) / x) < dxrel) { *!*/
        if (x < -0.5 && abs(x - (int) (x - 0.5) / x) < dxrel) {
          throw new ArithmeticException("Math Error: PRECISION");
        }

        /* The argument is so close to 0 that the result would overflow. */
        if (y < xsml) {
          if (x > 0) return Double.POSITIVE_INFINITY;
          return Double.NEGATIVE_INFINITY;
        }
        n = -n;
        for (i = 0; i < n; i++) {
          value /= (x + i);
        }
        return value;
      } else {
        /* gamma(x) for 2 <= x <= 10 */

        for (i = 1; i <= n; i++) {
          value *= (y + i);
        }
        return value;
      }
    } else {
      /* gamma(x) for	 y = |x| > 10. */
      if (x > xmax) {      /* Overflow */
        return Double.POSITIVE_INFINITY;
      }

      if (x < xmin) {      /* Underflow */
        return 0.;
      }

      if (y <= 50 && y == (int) y) { /* compute (n - 1)! */
        value = 1.;
        for (i = 2; i < y; i++) value *= i;
      } else { /* normal case */
        value = exp((y - 0.5) * log(y) - y + M_LN_SQRT_2PI + ((2 * y == 2 * y) ? stirlerr(y) : lgammacor(y)));
      }

      if (x > 0) return value;

      /*!* 	if (fabs((x - (int)(x - 0.5))/x) < dxrel){ *!*/
      if (abs((x - (int) (x - 0.5)) / x) < dxrel) {

        /* The answer is less than half precision because */
        /* the argument is too near a negative integer. */

        throw new ArithmeticException("Math Error: PRECISION");
      }

      /*!* 	sinpiy = sin(M_PI * y); *!*/
      double sinpiy = sin(M_PI * y);
      if (sinpiy == 0)    /* Negative integer arg - overflow */ return Double.POSITIVE_INFINITY;

      return -M_PI / (y * sinpiy * value);
    }
  }

  private static double stirlerr(double n) {
    final double S0 = 0.083333333333333333333, S1 = 0.00277777777777777777778, S2 = 0.00079365079365079365079365, S3 = 0.000595238095238095238095238, S4 = 0.0008417508417508417508417508;
    final double[] sferr_halves = { 0.0, /* n=0 - wrong, place holder only */
                                    0.1534264097200273452913848,  /* 0.5 */
                                    0.0810614667953272582196702,  /* 1.0 */
                                    0.0548141210519176538961390,  /* 1.5 */
                                    0.0413406959554092940938221,  /* 2.0 */
                                    0.03316287351993628748511048, /* 2.5 */
                                    0.02767792568499833914878929, /* 3.0 */
                                    0.02374616365629749597132920, /* 3.5 */
                                    0.02079067210376509311152277, /* 4.0 */
                                    0.01848845053267318523077934, /* 4.5 */
                                    0.01664469118982119216319487, /* 5.0 */
                                    0.01513497322191737887351255, /* 5.5 */
                                    0.01387612882307074799874573, /* 6.0 */
                                    0.01281046524292022692424986, /* 6.5 */
                                    0.01189670994589177009505572, /* 7.0 */
                                    0.01110455975820691732662991, /* 7.5 */
                                    0.010411265261972096497478567, /* 8.0 */
                                    0.009799416126158803298389475, /* 8.5 */
                                    0.009255462182712732917728637, /* 9.0 */
                                    0.008768700134139385462952823, /* 9.5 */
                                    0.008330563433362871256469318, /* 10.0 */
                                    0.007934114564314020547248100, /* 10.5 */
                                    0.007573675487951840794972024, /* 11.0 */
                                    0.007244554301320383179543912, /* 11.5 */
                                    0.006942840107209529865664152, /* 12.0 */
                                    0.006665247032707682442354394, /* 12.5 */
                                    0.006408994188004207068439631, /* 13.0 */
                                    0.006171712263039457647532867, /* 13.5 */
                                    0.005951370112758847735624416, /* 14.0 */
                                    0.005746216513010115682023589, /* 14.5 */
                                    0.005554733551962801371038690  /* 15.0 */ };
    double nn;
    if (n <= 15.0) {
      nn = n + n;
      if (nn == (int) nn) return (sferr_halves[(int) nn]);
      return (lgammafn(n + 1.) - (n + 0.5) * log(n) + n - M_LN_SQRT_2PI);
    }
    nn = n * n;
    if (n > 500) return ((S0 - S1 / nn) / n);
    if (n > 80) return ((S0 - (S1 - S2 / nn) / nn) / n);
    if (n > 35) return ((S0 - (S1 - (S2 - S3 / nn) / nn) / nn) / n);
    // 15 < n <= 35 :
    return ((S0 - (S1 - (S2 - (S3 - S4 / nn) / nn) / nn) / nn) / n);
  }

  private static double lgammafn(double x) {
    return lgammafn_sign(x, null);
  }

}
