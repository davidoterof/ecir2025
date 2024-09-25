package org.irlab.ecir25.metric;

import org.irlab.ecir25.util.Qrels;

public enum MetricEnum {
  MAP {
    @Override
    public Metric getNewInstance(final Qrels qrels) {
      return new MAP(qrels);
    }
  }, NDCG {
    @Override
    public Metric getNewInstance(final Qrels qrels) {
      return new NDCG(qrels);
    }
  };

  public abstract Metric getNewInstance(final Qrels qrels);
}
