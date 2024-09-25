package org.irlab.ecir25.nhst;

public interface CorrectionProcedure {
  double[] correct(double[] unadjustedPvalues, double alpha);
}
