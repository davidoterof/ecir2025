package org.irlab.ecir25.nhst;

import com.google.common.primitives.Doubles;
import org.jblas.DoubleMatrix;

import java.util.Arrays;
import java.util.Collections;
import java.util.concurrent.ThreadLocalRandom;

public final class RandomisedTukeyHSD extends MulticomparisonsTest {

  private final int permutations;

  public RandomisedTukeyHSD(int permutations) {
    this.permutations = permutations;
  }

  @Override
  protected double[] computePValues(DoubleMatrix data, double[] observedDifferences) {
    double[][] dataP = data.toArray2();
    double[] pvalues = new double[observedDifferences.length];
    double[] permutedDifferences = new double[permutations];
    for (int i = 0; i < permutations; i++) {
      double[][] permutedMatrix = getPermutation(dataP);
      double[] columMeans = computeColumMeans(permutedMatrix);
      double maxMean = Doubles.max(columMeans);
      double minMean = Doubles.min(columMeans);
      double permuttedDifference = maxMean - minMean;
      permutedDifferences[i] = permuttedDifference;
    }
    Arrays.sort(permutedDifferences);
    for (int i = 0; i < observedDifferences.length; i++) {
      for (int j = 0; j < permutations; j++) {
        if (permutedDifferences[j] >= Math.abs(observedDifferences[i])) {
          int count = permutations - j;
          pvalues[i] = (double) count / permutations;
          break;
        }
      }
    }

    return pvalues;
  }


  private double[] computeColumMeans(double[][] matrix) {
    int numberOfSystems = matrix[0].length;
    double[] means = new double[numberOfSystems];
    for (int col = 0; col < numberOfSystems; col++) {
      double sum = 0;
      for (double[] doubles : matrix) {
        sum += doubles[col];
      }
      means[col] = sum / matrix.length;
    }
    return means;
  }

  private double[][] getPermutation(double[][] matrix) {
    for (double[] row : matrix) {
      Collections.shuffle(Doubles.asList(row), ThreadLocalRandom.current());
    }
    return matrix;
  }
}
