package org.irlab.ecir25.models;

import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.irlab.ecir25.metric.Metric;
import org.irlab.ecir25.metric.MetricEnum;
import org.irlab.ecir25.util.ParallelArrays;
import org.irlab.ecir25.util.Qrels;
import org.irlab.ecir25.util.TRECParser;
import smile.classification.CustomLogitLogisticRegressor;

import java.nio.file.Path;
import java.util.*;
import java.util.Map.Entry;
import java.util.stream.Collectors;

public final class ProbabilityModel {
  public static final Random randomGenerator = new Random();
  private static final Logger LOG = LogManager.getLogger();
  private final Metric metric;
  private final Qrels qrels;
  private final Int2DoubleMap w1 = new Int2DoubleOpenHashMap();
  private final Int2DoubleMap w2 = new Int2DoubleOpenHashMap();
  private final Int2ObjectMap<double[][]> positions = new Int2ObjectOpenHashMap<>();
  private final Int2ObjectMap<int[]> relevantsInt = new Int2ObjectOpenHashMap<>();
  private final Int2ObjectMap<double[]> relevants = new Int2ObjectOpenHashMap<>();
  private final Int2ObjectMap<Int2ObjectMap<double[]>> probCache = new Int2ObjectOpenHashMap<>();
  private String runPath;

  public ProbabilityModel(Path pathToQrel, MetricEnum metric) {
    this.qrels = Qrels.fromFile(pathToQrel.toString());
    this.metric = Metric.getNewMetric(metric, this.qrels);
  }

  public void learnProbabilities(String pathTorun) {
    this.runPath = pathTorun;
    Map<Integer, List<String>> run = TRECParser.parseRun(pathTorun);
    for (Integer query : run.keySet()) {
      learn(query, run.get(query), qrels);
    }
  }

  public double[] sampleMetric(double change, int queries, int queriesImproved) {
    int[] sampledQueries = ArrayUtils.subarray(probCache.keySet().toIntArray(), 0, queries);
    return sampleMetric(change, queriesImproved, sampledQueries);
  }

  public double[] sampleMetric(double change, int queriesImproved, int[] queryIds) {
    Int2ObjectMap<int[]> sampledRelevants = new Int2ObjectOpenHashMap<>(queryIds.length);
    int count = 0;
    Set<Integer> improvedQueryIds = getRandomQueries(queryIds, queriesImproved);

    for (Integer queryId : queryIds) {
      if (improvedQueryIds.contains(queryId)) {
        sampledRelevants.put(queryId.intValue(), sampleRelevance(queryId, change));
      } else {
        sampledRelevants.put(queryId.intValue(), sampleRelevance(queryId, 0));
      }
      count++;
    }

    Map<Integer, Double> metricValues = metric.calculateFromRelevance(sampledRelevants);
    int[] keys = new int[metricValues.size()];
    double[] values = new double[metricValues.size()];
    count = 0;
    for (Entry<Integer, Double> entry : metricValues.entrySet()) {
      keys[count] = entry.getKey();
      values[count] = entry.getValue();
      count++;
    }
    ParallelArrays.sort(Arrays.stream(keys).asDoubleStream().toArray(), values, false);
    return values;
  }

  private Set<Integer> getRandomQueries(int[] queryIds, int number) {
    Set<Integer> selectedQueries = new HashSet<>();
    List<Integer> list = Arrays.stream(queryIds).boxed().collect(Collectors.toList());
    Collections.shuffle(list);
    for (int i = 0; i < number; i++) {
      selectedQueries.add(list.get(i));
    }
    return selectedQueries;
  }

  public void learn(int queryId, List<String> ranking, Qrels qrels) {
    double[] qRelevants = new double[ranking.size()];
    int[] qRelevantsInt = new int[ranking.size()];
    double[][] qPositions = new double[ranking.size()][1];
    probCache.put(queryId, new Int2ObjectOpenHashMap<>());
    int count = 0;
    for (String doc : ranking) {
      qRelevants[count] = qrels.getRelevance(queryId, doc);
      count++;
    }
    for (int i = 0; i < qPositions.length; i++) {
      qPositions[i][0] = i + 1d;
      qRelevantsInt[i] = (int) qRelevants[i];
    }

    relevants.put(queryId, qRelevants);
    relevantsInt.put(queryId, qRelevantsInt);
    positions.put(queryId, qPositions);

    CustomLogitLogisticRegressor regressor;
    try {
      regressor = CustomLogitLogisticRegressor.fit(qPositions, qRelevantsInt);
    } catch (IllegalArgumentException e) {
      // This happens when system retrieved no relevants for this topic.
      LOG.warn("Failed to fit regressor for query {} and run {}", queryId, runPath);
      return;
    }
    double[] w = regressor.getW();
    w1.put(queryId, w[0]);
    w2.put(queryId, w[1]);
  }

  public int[] sampleRelevance(int queryID, double change) {
    double[] tmp = new double[0];
    int cacheKey = (int) (change * 1000);
    synchronized (probCache) {
      Int2ObjectMap<double[]> cache = probCache.get(queryID);
      if (cache.get(cacheKey) == null) {
        boolean skip = false;
        double[] probabilities = new double[positions.get(queryID).length];
        double[] wchange = new double[2];
        double p = 1d + change;
        if (w1.get(queryID) > 0) {
          wchange[0] = p * w1.get(queryID);
        } else {
          wchange[0] = (1 / p) * w1.get(queryID);
        }
        if (w2.get(queryID) > 0) {
          wchange[1] = p * w2.get(queryID);
        } else {
          wchange[1] = (1 / p) * w2.get(queryID);
        }
        CustomLogitLogisticRegressor changedRegressor = null;
        try {
          changedRegressor = CustomLogitLogisticRegressor.dummyAlteredRegressor(positions.get(queryID),
                                                                                relevantsInt.get(queryID),
                                                                                wchange);
        } catch (IllegalArgumentException e) {
          // This happens when a ranking has no relevants. In that case,
          // the simulated systems will output non-relevants always.
          tmp = relevants.get(queryID);
          probCache.get(queryID).put(cacheKey, tmp);
          LOG.warn("Problems with regressor for query {} and run {}", queryID, runPath);
          skip = true; // We are not going to predict with the regressor for this query
        }
        if (!skip) {
          double[] posterior = new double[2];
          for (int i = 0; i < positions.get(queryID).length; i++) {
            changedRegressor.predict(new double[] { i + 1 }, posterior);
            probabilities[i] = posterior[1];
          }
          tmp = probabilities;
          probCache.get(queryID).put(cacheKey, probabilities);
        }
      } else {
        tmp = cache.get(cacheKey);
      }
    }

    try {
      tmp = probCache.get(queryID).get(cacheKey);
    } catch (java.lang.ArrayIndexOutOfBoundsException e) {
      LOG.error("Problems", e);
    }

    int[] sample = new int[tmp.length];
    for (int i = 0; i < sample.length; i++) {
      if (randomGenerator.nextDouble() <= tmp[i]) {
        sample[i] = 1;
      } else {
        sample[i] = 0;
      }
    }
    return sample;
  }
}
