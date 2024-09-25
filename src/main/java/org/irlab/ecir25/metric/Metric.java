package org.irlab.ecir25.metric;

import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import org.irlab.ecir25.util.Qrels;

import java.util.List;
import java.util.Map;

public interface Metric {

  static Metric getNewMetric(final MetricEnum metric, final Qrels qrels) {
    return metric.getNewInstance(qrels);
  }

  Map<Integer, Double> calculateFromRelevance(Int2ObjectMap<int[]> rankingsRelevance);

  double calculatePerQueryFromRelevance(int[] rankingRelevance, int queryId);

  double computeForTopic(final int topicID, List<String> docsRanking);
}
