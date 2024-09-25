package org.irlab.ecir25.nhst;

import org.irlab.ecir25.util.ParallelArrays;

import java.util.Arrays;
import java.util.stream.IntStream;

public class HolmCorrection implements CorrectionProcedure {
  @Override
  public double[] correct(double[] unadjustedPvalues, double alpha) {
    int total_comparisons = unadjustedPvalues.length;
    double[] unadjusted = Arrays.copyOf(unadjustedPvalues, unadjustedPvalues.length);
    double[] positions = IntStream.range(0, unadjusted.length).mapToDouble(val -> (double) val).toArray();
    ParallelArrays.sort(unadjusted, positions, true);
    double[] newPvalues = new double[total_comparisons];
    for (int i = 0; i < total_comparisons; i++) {
      double pvalue = unadjusted[i];
      if (pvalue <= (alpha / (total_comparisons - i))) {
        newPvalues[(int) positions[i]] = 0;
      } else {
        newPvalues[(int) positions[i]] = 1;
      }
    }
    return newPvalues;
  }
}

