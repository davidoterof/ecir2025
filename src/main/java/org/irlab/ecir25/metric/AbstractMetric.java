package org.irlab.ecir25.metric;

import com.koloboke.collect.map.hash.HashIntDoubleMaps;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import org.irlab.ecir25.util.Qrels;

import java.util.Map;

abstract class AbstractMetric implements Metric {

  protected final Qrels qrels;

  protected AbstractMetric(Qrels qrels) {
    this.qrels = qrels;
  }

  protected abstract String getName();

  @Override
  public final String toString() {
    return getName();
  }

  @Override
  public Map<Integer, Double> calculateFromRelevance(Int2ObjectMap<int[]> rankingsRelevance) {
    final Map<Integer, Double> maps = HashIntDoubleMaps.newUpdatableMap(rankingsRelevance.size());
    rankingsRelevance.forEach((queryId, relevanceValues) -> {
      maps.put(queryId, calculatePerQueryFromRelevance(relevanceValues, queryId));
    });

    return maps;
  }
}
