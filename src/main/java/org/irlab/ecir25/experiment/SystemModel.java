package org.irlab.ecir25.experiment;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.irlab.ecir25.models.ProbabilityModel;
import org.irlab.ecir25.nhst.*;
import org.jblas.DoubleMatrix;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.*;

public class SystemModel {
  private final ProbabilityModel model;
  private final int numberOfSimulationsPerSystem;
  private final int threads;
  private final String runPath;

  public SystemModel(int numberOfSimulationsPerSystem, int threads, ProbabilityModel model, String runPath) {
    this.numberOfSimulationsPerSystem = numberOfSimulationsPerSystem;
    this.threads = threads;
    this.model = model;
    this.runPath = runPath;
  }

  public void fitSystem() {
    model.learnProbabilities(runPath);
  }

  public double[][] compareUnderProbabilisticModel(int queries, int numberOfSystems, double es, double alpha,
                                                   int number_of_tests) {
    try {
      double[][] results = new double[number_of_tests][numberOfSimulationsPerSystem];
      ExecutorService pool = Executors.newFixedThreadPool(threads);

      List<Future<Pair<Integer, int[]>>> futures = new ArrayList<>();
      for (int i = 0; i < numberOfSimulationsPerSystem; i++) {
        futures.add(pool.submit(new CompareTestsThread(model, i, queries, numberOfSystems, es, alpha)));
      }

      for (Future<Pair<Integer, int[]>> future : futures) {
        int j = future.get().getLeft();
        int[] result = future.get().getRight();
        for (int i = 0; i < number_of_tests; i++) {
          results[i][j] = result[i];
        }
      }

      pool.shutdown();
      pool.awaitTermination(5, TimeUnit.MINUTES);
      pool.close();
      return results;
    } catch (ExecutionException | InterruptedException e) {
      throw new RuntimeException(e);
    }
  }

  private static class CompareTestsThread implements Callable<Pair<Integer, int[]>> {
    private final Logger LOG = LogManager.getLogger();

    private final ProbabilityModel model;
    private final int i;
    private final int queries;
    private final int numberOfSystems;
    private final double es;
    private final double alpha;

    public CompareTestsThread(ProbabilityModel model, int i, int queries, int numberOfSystems, double es,
                              double alpha) {
      this.model = model;
      this.i = i;
      this.queries = queries;
      this.numberOfSystems = numberOfSystems;
      this.es = es;
      this.alpha = alpha;
    }

    @Override
    public Pair<Integer, int[]> call() {
      try {

        int posFitted = (numberOfSystems / 2) + 1;
        int limit = numberOfSystems - posFitted;
        int upperLimit;
        if (numberOfSystems % 2 == 0) {
          upperLimit = limit + 1;
        } else {
          upperLimit = limit;
        }

        DoubleMatrix scoreMatrix = null;
        for (int i = -limit; i <= upperLimit; i++) {
          double[] scores = model.sampleMetric(i * es, queries, queries);
          if (scoreMatrix == null) {
            scoreMatrix = new DoubleMatrix(scores);
          } else {
            scoreMatrix = DoubleMatrix.concatHorizontally(scoreMatrix, new DoubleMatrix(scores));
          }
        }

        double[] unadjustedTTestPValues = new UnadjustedTTest().test(scoreMatrix);
        int unadjustedTTestRejected = 0;
        for (double pvalue : unadjustedTTestPValues) {
          if (pvalue <= alpha) unadjustedTTestRejected++;
        }
        int unadjustedTTestError = computeMetric(unadjustedTTestPValues.length, unadjustedTTestRejected, es);

        double[] bonferroniTTestPValues = new BonferroniCorrection().correct(unadjustedTTestPValues, alpha);
        int bonferroniTTestRejected = 0;
        for (double pvalue : bonferroniTTestPValues) {
          if (pvalue <= alpha) bonferroniTTestRejected++;
        }
        int bonferroniTTestError = computeMetric(bonferroniTTestPValues.length, bonferroniTTestRejected, es);

        double[] holmTTestPValues = new HolmCorrection().correct(unadjustedTTestPValues, alpha);
        int holmTTestRejected = 0;
        for (double pvalue : holmTTestPValues) {
          if (pvalue <= alpha) holmTTestRejected++;
        }
        int holmTTestError = computeMetric(holmTTestPValues.length, holmTTestRejected, es);

        double[] bhTTestPValues = new BHCorrection().correct(unadjustedTTestPValues, alpha);
        int bhTTestRejected = 0;
        for (double pvalue : bhTTestPValues) {
          if (pvalue <= alpha) bhTTestRejected++;
        }
        int bhTTestError = computeMetric(bhTTestPValues.length, bhTTestRejected, es);

        double[] byTTestPvalues = new BYCorrection().correct(unadjustedTTestPValues, alpha);
        int byTTestRejected = 0;
        for (double pvalue : byTTestPvalues) {
          if (pvalue <= alpha) byTTestRejected++;
        }
        int byTTestError = computeMetric(byTTestPvalues.length, byTTestRejected, es);

        double[] tukeyPvalues = new RandomisedTukeyHSD(100_000).test(scoreMatrix);
        int tukeyRejected = 0;
        for (double pvalue : tukeyPvalues) {
          if (pvalue <= alpha) tukeyRejected++;
        }
        int tukeyError = computeMetric(tukeyPvalues.length, tukeyRejected, es);

        double[] unadjustedWilcoxonPValues = new UnadjustedWilcoxon().test(scoreMatrix);
        int unadjustedWilcoxonRejected = 0;
        for (double pvalue : unadjustedWilcoxonPValues) {
          if (pvalue <= alpha) unadjustedWilcoxonRejected++;
        }
        int unadjustedWilcoxonError = computeMetric(unadjustedWilcoxonPValues.length, unadjustedWilcoxonRejected, es);

        double[] boferroniWilcoxonPValues = new BonferroniCorrection().correct(unadjustedWilcoxonPValues, alpha);
        int bonferroniWilcoxonRejected = 0;
        for (double pvalue : boferroniWilcoxonPValues) {
          if (pvalue <= alpha) bonferroniWilcoxonRejected++;
        }
        int bonferroniWilcoxonError = computeMetric(boferroniWilcoxonPValues.length, bonferroniWilcoxonRejected, es);

        double[] holmWilcoxonPValues = new HolmCorrection().correct(unadjustedWilcoxonPValues, alpha);
        int holmWilcoxonRejected = 0;
        for (double pvalue : holmWilcoxonPValues) {
          if (pvalue <= alpha) holmWilcoxonRejected++;
        }
        int holmWilcoxonError = computeMetric(holmWilcoxonPValues.length, holmWilcoxonRejected, es);

        double[] bhWilcoxonPValues = new BHCorrection().correct(unadjustedWilcoxonPValues, alpha);
        int bhWilcoxonRejected = 0;
        for (double pvalue : bhWilcoxonPValues) {
          if (pvalue <= alpha) bhWilcoxonRejected++;
        }
        int bhWilcoxonError = computeMetric(bhWilcoxonPValues.length, bhWilcoxonRejected, es);

        double[] byWilcoxonPValues = new BYCorrection().correct(unadjustedWilcoxonPValues, alpha);
        int byWilcoxonRejected = 0;
        for (double pvalue : byWilcoxonPValues) {
          if (pvalue <= alpha) byWilcoxonRejected++;
        }
        int byWilcoxonError = computeMetric(byWilcoxonPValues.length, byWilcoxonRejected, es);

        double[] anovapvalues = new TwoWayANOVA().test(scoreMatrix);
        int anovaRejected = 0;
        for (double pvalue : anovapvalues) {
          if (pvalue <= alpha) anovaRejected++;
        }
        int anovaError = computeMetric(anovapvalues.length, anovaRejected, es);

        return Pair.of(i,
                       new int[] { unadjustedTTestError,
                                   bonferroniTTestError,
                                   holmTTestError,
                                   bhTTestError,
                                   byTTestError,
                                   tukeyError,
                                   unadjustedWilcoxonError,
                                   bonferroniWilcoxonError,
                                   holmWilcoxonError,
                                   bhWilcoxonError,
                                   byWilcoxonError,
                                   anovaError });
      } catch (Exception e) {
        LOG.warn("Failed to simulate comparison, including all p-values 1", e);
        return Pair.of(i, new int[] { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 });
      }
    }

    private int computeMetric(int numHypotheses, int rejected, double es) {
      if (es == 0) {
        return rejected > 0 ? 1 : 0;
      } else {
        return numHypotheses == rejected ? 1 : 0;
      }
    }
  }
}
