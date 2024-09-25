package org.irlab.ecir25.nhst;

import org.apache.commons.math3.stat.StatUtils;
import org.jblas.DoubleMatrix;

// Implementation of ANOVA + TukeyHSD correction.
public final class TwoWayANOVA extends MulticomparisonsTest {

  @Override
  protected double[] computePValues(DoubleMatrix data, double[] observedDifferences) {
    int number_of_systems = data.columns;
    int number_of_topics = data.rows;

    int dfe = (number_of_topics * number_of_systems) - number_of_systems;

    double ss_error = 0;

    for (int i = 0; i < data.columns; i++) {
      double[] column = data.getColumn(i).data;
      ss_error += StatUtils.variance(column) * (number_of_systems - 1);
    }

    double ms_error = ss_error / dfe;

    double[] pvalues = new double[observedDifferences.length];

    for (int i = 0; i < observedDifferences.length; i++) {
      double observedDifference = observedDifferences[i];

      double stand_err = Math.sqrt(ms_error / number_of_topics);

      double tukeyhsd_statistic = Math.abs(observedDifference) / stand_err;

      double pvalue = Tukey.cumulative(tukeyhsd_statistic, 1, number_of_systems, dfe, false, false);
      pvalues[i] = pvalue;
    }

    return pvalues;
  }
}
