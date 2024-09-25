package org.irlab.ecir25.nhst;

public class BonferroniCorrection implements CorrectionProcedure {
  @Override
  public double[] correct(double[] unadjustedPvalues, double alpha) {
    int total_comparisons = unadjustedPvalues.length;
    double[] correctedPValues = new double[total_comparisons];
    for (int i = 0; i < total_comparisons; i++) {
      if (unadjustedPvalues[i] <= (alpha / total_comparisons)) {
        correctedPValues[i] = 0;
      } else {
        correctedPValues[i] = 1;
      }
    }
    return correctedPValues;
  }
}
