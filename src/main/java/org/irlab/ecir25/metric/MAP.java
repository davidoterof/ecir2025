package org.irlab.ecir25.metric;

import org.irlab.ecir25.util.Qrels;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.List;

final class MAP extends AbstractMetric {

  MAP(final Qrels qrels) {
    super(qrels);
  }

  @Override
  public double calculatePerQueryFromRelevance(int[] rankingRelevance, final int topicID) {
    double ap = 0;
    double relevants = 0;
    for (int i = 0; i < rankingRelevance.length; i++) {
      if (rankingRelevance[i] > 0) {
        relevants++;
        double p = relevants / (i + 1);
        ap += p;
      }
    }

    final long totalTopicRelevants = qrels.totalRelevantsForTopic(topicID);
    ap = ap == 0 ? 0 : ap / (double) totalTopicRelevants;
    return Math.min(ap, 1);
  }

  @Override
  public double computeForTopic(final int topicID, final List<String> docsRanking) {
    double ap = 0;
    double relevants = 0;
    double count = 0;
    for (final String doc : docsRanking) {
      count++;
      if (qrels.isDocumentRelevant(topicID, doc)) {
        relevants++;
        double p = relevants / count;
        ap = ap + p;
      }
    }

    final long totalTopicRelevants = qrels.totalRelevantsForTopic(topicID);
    ap = ap == 0 ? 0 : ap / (double) totalTopicRelevants;
    return BigDecimal.valueOf(ap).setScale(5, RoundingMode.HALF_UP).doubleValue();
  }


  @Override
  public String getName() {
    return "MAP";
  }
}
