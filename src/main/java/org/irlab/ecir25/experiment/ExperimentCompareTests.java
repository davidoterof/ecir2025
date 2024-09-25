package org.irlab.ecir25.experiment;

import org.apache.commons.cli.CommandLine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.irlab.ecir25.metric.MetricEnum;
import org.irlab.ecir25.models.ProbabilityModel;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.nio.file.StandardOpenOption;

public class ExperimentCompareTests {

  private static final Logger LOG = LogManager.getLogger();

  private static final int SIMULATIONS = 1000;

  public static void run(CommandLine args) {
    int THREADS = Runtime.getRuntime().availableProcessors();
    LOG.info("Running ExperimentCompareTests with {} threads", THREADS);
    String collection = args.getOptionValue("c");
    Path runsPath = Path.of(args.getOptionValue("r"));
    Path qrelsPath = Path.of(args.getOptionValue("q"));
    Path outputFolderPath = Path.of(args.getOptionValue("o"));

    MetricEnum[] metrics = new MetricEnum[] { MetricEnum.MAP, MetricEnum.NDCG };
    double[] alphas = new double[] { 0.05 };
    int[] samplesizes = new int[] { 10, 30, 50 };
    int[] numbersOfSystems = new int[] { 3, 5, 10 };

    for (MetricEnum metric : metrics) {
      for (double alpha : alphas) {
        for (int numberOfSystems : numbersOfSystems) {
          for (int sampleSize : samplesizes) {
            String folder = String.format("logreg_%s_alpha-%1.2f_%s_systems-%d_topics-%d_sims_%d",
                                          collection,
                                          alpha,
                                          metric,
                                          numberOfSystems,
                                          sampleSize,
                                          SIMULATIONS);
            Path experimentOutput = outputFolderPath.resolve(folder);
            experimentOutput.toFile().mkdir();

            String[] files = runsPath.toFile().list();
            for (String runName : files) {
              LOG.info("PROCESSING {} from {} with metric {}, sample size {}, and number of systems {}",
                       runName,
                       collection,
                       metric,
                       sampleSize,
                       numberOfSystems);
              String runPath = runsPath + File.separator + runName;

              SystemModel systemModel = new SystemModel(SIMULATIONS,
                                                        THREADS,
                                                        new ProbabilityModel(qrelsPath, metric),
                                                        runPath);
              systemModel.fitSystem();

              for (double es = 0.00d; es <= 0.255; es += 0.005d) {
                boolean created = createResultsFile(experimentOutput, runName, es);

                if (!created) {
                  LOG.warn("SKIPPING {} from TREC-{} with metric {}, sample size {}, number of systems {} and es {}",
                           runName,
                           collection,
                           metric,
                           sampleSize,
                           numberOfSystems,
                           es);
                  continue;
                }

                LOG.info("COMPUTING {} from TREC-{} with metric {}, sample size {}, number of systems {} and es {}",
                         runName,
                         collection,
                         metric,
                         sampleSize,
                         numberOfSystems,
                         es);

                int number_of_tests = 12;

                double[][] results = systemModel.compareUnderProbabilisticModel(sampleSize,
                                                                                numberOfSystems,
                                                                                es,
                                                                                alpha,
                                                                                number_of_tests);

                try (final PrintWriter writer = new PrintWriter(Files.newBufferedWriter(getTempOutputFile(
                    experimentOutput,
                    runName,
                    es), StandardOpenOption.CREATE, StandardOpenOption.APPEND))) {
                  writer.println(
                      "UnadjustedTTest\tBonferroniTTest\tHolmTTest\tBHTTest\tBYTTest\tTukey\tUnadjustedWilcoxon\tBonferroniWilcoxon\tHolmWilcoxon\tBHWilcoxon\tBYWilcoxon\tANOVA");
                  double[] power = new double[number_of_tests];

                  for (int j = 0; j < SIMULATIONS; j++) {
                    for (int k = 0; k < number_of_tests; k++) {
                      power[k] += results[k][j];
                    }
                  }

                  for (int i = 0; i < number_of_tests - 1; i++) {
                    writer.print(power[i] / SIMULATIONS);
                    writer.print("\t");
                  }
                  writer.print(power[number_of_tests - 1] / SIMULATIONS);
                  writer.println();
                  writer.flush();

                  Path outputFile = getOutputFile(experimentOutput, runName, es);
                  Path tempFile = getTempOutputFile(experimentOutput, runName, es);
                  Files.move(tempFile, outputFile, StandardCopyOption.ATOMIC_MOVE);
                } catch (IOException e) {
                  LOG.error("problem in experiments", e);
                }
              }
            }
          }
        }
      }
    }
  }

  private static boolean createResultsFile(Path outputFolder, String runFile, double es) {
    final Path outputPath = getOutputFile(outputFolder, runFile, es);
    final Path tempPath = getTempOutputFile(outputFolder, runFile, es);

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

  private static String getOutputFilename(String runFile, double es) {
    return String.format("%s-es-%.3f.tsv", runFile, es);
  }

  private static Path getOutputFile(Path outputFolder, String runFile, double es) {
    final String filename = getOutputFilename(runFile, es);
    return outputFolder.resolve(filename);
  }

  private static Path getTempOutputFile(Path outputFolder, String runFile, double es) {
    final String filename = getOutputFilename(runFile, es);
    return outputFolder.resolve(filename.replace(".tsv", ".tmp"));
  }
}
