package org.irlab.ecir25.nhst;

import org.jblas.DoubleMatrix;

public abstract class MulticomparisonsTest {

  public double[] test(DoubleMatrix data) {
    double[] observedDifferences = computeMeansDifferences(data);
    return computePValues(data, observedDifferences);
  }

  protected abstract double[] computePValues(DoubleMatrix data, double[] observedDifferences);

  protected double[] computeMeansDifferences(DoubleMatrix data) {
    double[] means = data.columnMeans().data;
    int systems = means.length;
    int totalComparisons = ((systems - 1) * systems) / 2;
    double[] meansDifferences = new double[totalComparisons];
    int count = 0;
    for (int i = 0; i < systems; i++) {
      for (int j = i + 1; j < systems; j++) {
        double difference = means[j] - means[i];
        meansDifferences[count] = difference;
        count++;
      }
    }

    return meansDifferences;
  }
}
