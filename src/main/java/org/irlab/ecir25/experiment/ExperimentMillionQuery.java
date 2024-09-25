package org.irlab.ecir25.experiment;

import org.apache.commons.cli.CommandLine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.irlab.ecir25.metric.Metric;
import org.irlab.ecir25.metric.MetricEnum;
import org.irlab.ecir25.nhst.*;
import org.irlab.ecir25.util.Qrels;
import org.irlab.ecir25.util.Runs;
import org.jblas.DoubleMatrix;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.nio.file.StandardOpenOption;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;

public class ExperimentMillionQuery {

  private static final Logger LOG = LogManager.getLogger();

  private final static int MIN_SAMPLES = 1;
  private final static int MAX_SAMPLES = 2000;
  private final static int MIN_TOPICS = 10;
  private final static int MAX_TOPICS = 100;

  public static void run(CommandLine args) {
    Path outputFolderPath = Path.of(args.getOptionValue("o"));
    Path runsPath = Path.of(args.getOptionValue("r"));
    String qrelsPath = args.getOptionValue("q");

    LOG.info("Reading runs from disk");
    File[] runsFiles = runsPath.toFile().listFiles();
    List<String> runsPaths = Arrays.stream(runsFiles)
                                   .map(File::getAbsolutePath)
                                   .sorted(String.CASE_INSENSITIVE_ORDER)
                                   .collect(Collectors.toList());
    Runs runs = Runs.fromListOfPaths(runsPaths);
    LOG.info("Runs read");

    Qrels qrels = Qrels.fromFile(qrelsPath);
    Metric refMetric = Metric.getNewMetric(MetricEnum.MAP, qrels);

    // Matrix of scores: rows are topics, columns are systems,
    DoubleMatrix scoreMatrix = getMatrix(refMetric, runs, qrels);

    double[] meanScores = scoreMatrix.columnMeans().data;

    double[] goldPValues = new double[(meanScores.length * (meanScores.length - 1)) / 2];

    int count = 0;
    for (int i = 0; i < meanScores.length; i++) {
      for (int j = i + 1; j < meanScores.length; j++) {
        double s1 = Math.min(meanScores[i], meanScores[j]);
        double s2 = Math.max(meanScores[i], meanScores[j]);
        double diff = ((s2 - s1) * 100) / s1;
        if (diff >= 0.5d) {
          goldPValues[count] = 0;
        } else {
          goldPValues[count] = 1;
        }
        count++;
      }
    }

    // Write gold pvalues to disk.
    try (final PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outputFolderPath.resolve("gold.tsv"),
                                                                            StandardOpenOption.CREATE))) {
      writer.print("gold\t");
      for (double pvalue : goldPValues) {
        writer.print(pvalue);
        writer.print("\t");
      }
      writer.println();
    } catch (IOException e) {
      throw new RuntimeException(e);
    }

    ExecutorService pool = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
    for (int numberOfQueries = MIN_TOPICS; numberOfQueries <= MAX_TOPICS; numberOfQueries++) {
      for (int sample = MIN_SAMPLES; sample <= MAX_SAMPLES; sample++) {
        // Shuffle and sample topics.
        List<Integer> topics = new ArrayList<>(qrels.getTopics());
        Collections.shuffle(topics);
        List<Integer> sampledTopics = topics.subList(0, numberOfQueries);

        // Create metric with sampled topics.
        Qrels newQrels = qrels.filterTopics(sampledTopics);
        Metric testMetric = Metric.getNewMetric(MetricEnum.MAP, newQrels);

        // Obtain new scores.
        DoubleMatrix testMatrix = getMatrix(testMetric, runs, newQrels);

        pool.submit(new ExperimentThread(numberOfQueries, sample, outputFolderPath, testMatrix));
      }
    }
    pool.shutdown();
    try {
      pool.awaitTermination(5, TimeUnit.MINUTES);
    } catch (InterruptedException e) {
      throw new RuntimeException(e);
    }
  }

  private static DoubleMatrix getMatrix(Metric metric, Runs runs, Qrels qrels) {
    List<Integer> queries = qrels.getTopics();
    DoubleMatrix scoreMatrix = null;
    for (Map<Integer, List<String>> system : runs.getRuns()) {
      double[] runScores = evaluateSystem(system, queries, metric);
      scoreMatrix = appendToScoreMatrix(scoreMatrix, runScores);
    }
    return scoreMatrix;
  }

  private static double[] evaluateSystem(Map<Integer, List<String>> system, List<Integer> queries, Metric metric) {
    double[] runScores = new double[queries.size()];
    int j = 0;
    for (int queryID : queries) {
      List<String> docs = system.get(queryID);
      runScores[j] = evaluateSystemForQuery(docs, queryID, metric);
      j++;
    }
    return runScores;
  }

  private static double evaluateSystemForQuery(List<String> docs, int queryID, Metric metric) {
    return (docs == null || docs.isEmpty()) ? 0 : metric.computeForTopic(queryID, docs);
  }

  private static DoubleMatrix appendToScoreMatrix(DoubleMatrix scoreMatrix, double[] runScores) {
    if (scoreMatrix == null) {
      return new DoubleMatrix(runScores);
    } else {
      return DoubleMatrix.concatHorizontally(scoreMatrix, new DoubleMatrix(runScores));
    }
  }

  private static class ExperimentThread implements Runnable {

    private final int numberOfQueries;
    private final int sample;
    private final Path outputFolderPath;
    private final DoubleMatrix testMatrix;

    private ExperimentThread(int numberOfQueries, int sample, Path outputFolderPath, DoubleMatrix testMatrix) {
      this.numberOfQueries = numberOfQueries;
      this.sample = sample;
      this.outputFolderPath = outputFolderPath;
      this.testMatrix = testMatrix;
    }

    private String getOutputFilename() {
      return String.format("num-queries-%d-sample-%d.tsv", numberOfQueries, sample);
    }

    private Path getOutputFile() {
      return outputFolderPath.resolve(getOutputFilename());
    }

    private Path getTempOutputFile() {
      return outputFolderPath.resolve(getOutputFilename().replace(".tsv", ".tmp"));
    }

    private boolean createOutputFile() {
      final Path outputPath = getOutputFile();
      final Path tempPath = getTempOutputFile();

      if (Files.exists(outputPath) || Files.exists(tempPath)) {
        return false;
      }

      try {
        Files.createFile(tempPath);
      } catch (IOException e) {
        return false;
      }

      return true;
    }

    @Override
    public void run() {
      boolean created = createOutputFile();
      if (!created) {
        LOG.warn("Ignoring num-queries-{}-sample-{}", numberOfQueries, sample);
        return;
      }

      LOG.info("Computing num-queries-{}-sample-{}", numberOfQueries, sample);

      double[] tTestPValues = new UnadjustedTTest().test(testMatrix);
      double[] tTestBonferroniPValues = new BonferroniCorrection().correct(tTestPValues, 0.05d);
      double[] tTestHolmPValues = new HolmCorrection().correct(tTestPValues, 0.05d);
      double[] tTestBHPValues = new BHCorrection().correct(tTestPValues, 0.05d);
      double[] tTestBYPValues = new BYCorrection().correct(tTestPValues, 0.05d);

      double[] tukeyPValues = new RandomisedTukeyHSD(100_000).test(testMatrix);

      double[] wilcoxonPValues = new UnadjustedWilcoxon().test(testMatrix);
      double[] wilcoxonBonferroniPValues = new BonferroniCorrection().correct(wilcoxonPValues, 0.05d);
      double[] wilcoxonHolmPValues = new HolmCorrection().correct(wilcoxonPValues, 0.05d);
      double[] wilcoxonBHPValues = new BHCorrection().correct(wilcoxonPValues, 0.05d);
      double[] wilcoxonBYPValues = new BYCorrection().correct(wilcoxonPValues, 0.05d);

      try (final PrintWriter writer = new PrintWriter(Files.newBufferedWriter(getTempOutputFile(),
                                                                              StandardOpenOption.CREATE,
                                                                              StandardOpenOption.APPEND))) {
        writer.print("unadjusted t-test\t");
        for (double pvalue : tTestPValues) {
          writer.print(pvalue);
          writer.print("\t");
        }
        writer.println();

        writer.print("t-test bonferroni\t");
        for (double pvalue : tTestBonferroniPValues) {
          writer.print(pvalue);
          writer.print("\t");
        }
        writer.println();

        writer.print("t-test holm\t");
        for (double pvalue : tTestHolmPValues) {
          writer.print(pvalue);
          writer.print("\t");
        }
        writer.println();

        writer.print("t-test bh\t");
        for (double pvalue : tTestBHPValues) {
          writer.print(pvalue);
          writer.print("\t");
        }
        writer.println();

        writer.print("t-test by\t");
        for (double pvalue : tTestBYPValues) {
          writer.print(pvalue);
          writer.print("\t");
        }
        writer.println();

        writer.print("unadjusted wilcoxon\t");
        for (double pvalue : wilcoxonPValues) {
          writer.print(pvalue);
          writer.print("\t");
        }
        writer.println();

        writer.print("wilcoxon bonferroni\t");
        for (double pvalue : wilcoxonBonferroniPValues) {
          writer.print(pvalue);
          writer.print("\t");
        }
        writer.println();

        writer.print("wilcoxon holm\t");
        for (double pvalue : wilcoxonHolmPValues) {
          writer.print(pvalue);
          writer.print("\t");
        }
        writer.println();

        writer.print("wilcoxon bh\t");
        for (double pvalue : wilcoxonBHPValues) {
          writer.print(pvalue);
          writer.print("\t");
        }
        writer.println();

        writer.print("wilcoxon by\t");
        for (double pvalue : wilcoxonBYPValues) {
          writer.print(pvalue);
          writer.print("\t");
        }
        writer.println();

        writer.print("tukey\t");
        for (double pvalue : tukeyPValues) {
          writer.print(pvalue);
          writer.print("\t");
        }
        writer.println();

      } catch (IOException e) {
        throw new RuntimeException(e);
      }

      try {
        Files.move(getTempOutputFile(), getOutputFile(), StandardCopyOption.ATOMIC_MOVE);
      } catch (IOException e) {
        throw new RuntimeException(e);
      }
    }
  }
}
