package org.irlab.ecir25.metric;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.util.FastMath;
import org.irlab.ecir25.util.Qrels;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.Arrays;
import java.util.List;

final class NDCG extends AbstractMetric {

  NDCG(Qrels qrels) {
    super(qrels);
  }

  private double computeIDCG(final int queryID) {
    double[] values = qrels.getForTopic(queryID).values().stream().mapToDouble(Double::doubleValue).toArray();
    Arrays.sort(values);
    ArrayUtils.reverse(values);

    double idcg = 0d;
    for (int i = 0; i < values.length; i++) {
      final double rel = values[i];
      idcg = idcg + rel / (Math.log(i + 2) / Math.log(2));
    }

    return BigDecimal.valueOf(idcg).setScale(5, RoundingMode.HALF_UP).doubleValue();
  }

  @Override
  public double calculatePerQueryFromRelevance(int[] rankingRelevance, int queryId) {
    double idcg = computeIDCG(queryId);
    if (idcg == 0) {
      return 0;
    }
    double dcg = 0;
    for (int i = 0; i < rankingRelevance.length; i++) {
      double gain = rankingRelevance[i] / FastMath.log(2, i + 2);
      dcg += gain;
    }
    return Math.min(dcg / idcg, 1);
  }

  @Override
  public double computeForTopic(final int topicID, final List<String> docsRanking) {
    double idcg = computeIDCG(topicID);
    if (idcg == 0d) {
      return 0d;
    }
    double dcg = computeDCG(docsRanking, topicID);
    return Math.min(dcg / idcg, 1);
  }

  private double computeDCG(final List<String> docsRanking, final int queryID) {
    double dcg = 0;
    for (int i = 0; i < docsRanking.size(); i++) {
      String doc = docsRanking.get(i);
      double relevance = qrels.getRelevance(queryID, doc);
      double gain = relevance / FastMath.log(2, i + 2);
      dcg = dcg + gain;
    }
    return BigDecimal.valueOf(dcg).setScale(5, RoundingMode.HALF_UP).doubleValue();
  }

  @Override
  protected String getName() {
    return "NDCG";
  }
}
