package org.irlab.ecir25.nhst;

import org.apache.commons.math3.stat.inference.TTest;
import org.jblas.DoubleMatrix;

public final class UnadjustedTTest extends MulticomparisonsTest {

  @Override
  protected double[] computePValues(DoubleMatrix data, double[] observedDifferences) {
    int systems = data.columns;
    int totalComparisons = ((systems - 1) * systems) / 2;
    double[] pvalues = new double[totalComparisons];
    int count = 0;
    for (int i = 0; i < systems; i++) {
      for (int j = i + 1; j < systems; j++) {
        double[] sys1Scores = data.getColumn(i).data;
        double[] sys2Scores = data.getColumn(j).data;
        double pvalue = new TTest().pairedTTest(sys1Scores, sys2Scores);
        if (Double.isNaN(pvalue)) {
          pvalues[count] = 1d;
        } else {
          pvalues[count] = pvalue;
        }

        count++;
      }
    }
    return pvalues;
  }
}
