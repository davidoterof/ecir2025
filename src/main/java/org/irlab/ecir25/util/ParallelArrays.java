package org.irlab.ecir25.util;

import java.util.AbstractList;
import java.util.Comparator;

public final class ParallelArrays extends AbstractList<double[]> {

  private final double[] array1;
  private final double[] array2;

  public ParallelArrays(double[] array1, double[] array2) {
    if (array1.length != array2.length) throw new IllegalArgumentException();
    this.array1 = array1;
    this.array2 = array2;
  }

  public static void sort(double[] key, double[] array, boolean ascending) {
    if (ascending) {
      new ParallelArrays(key, array).sort(Comparator.comparingDouble(arr -> ((double[]) arr)[0]));
    } else {
      new ParallelArrays(key, array).sort(Comparator.comparingDouble(arr -> ((double[]) arr)[0]).reversed());
    }
  }

  public double[] get(int i) {
    return new double[] { array1[i], array2[i] };
  }

  public int size() {
    return array1.length;
  }

  public double[] set(int i, double[] a) {
    if (a.length != 2) throw new IllegalArgumentException();
    double[] b = get(i);
    array1[i] = a[0];
    array2[i] = a[1];
    return b;
  }

}
