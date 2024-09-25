package org.irlab.ecir25.util;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public final class Qrels {

  private final Map<Integer, Map<String, Double>> qrels;

  private final String name;

  private Qrels(final Map<Integer, Map<String, Double>> qrels) {
    this.qrels = qrels;
    this.name = "";
  }

  public static Qrels fromFile(final String path) {
    return new Qrels(TRECParser.parseQrels(path));
  }

  public Qrels filterTopics(final List<Integer> topics) {
    return new Qrels(qrels.entrySet()
                          .stream()
                          .filter(e -> topics.contains(e.getKey()))
                          .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue)));
  }

  public double getRelevance(final int topicID, final String document) {
    if (!qrels.containsKey(topicID)) {
      return 0d;
    }

    if (!qrels.get(topicID).containsKey(document)) {
      return 0d;
    }

    return qrels.get(topicID).get(document);
  }

  public List<Integer> getTopics() {
    return List.copyOf(qrels.keySet());
  }

  public Map<String, Double> getForTopic(final int topicID) {
    return qrels.get(topicID);
  }

  public int totalRelevantsForTopic(final int topicID) {
    return (int) getForTopic(topicID).entrySet().stream().filter(e -> !e.getValue().equals(0d)).count();
  }

  public boolean isDocumentRelevant(final int topicID, final String document) {
    return getRelevance(topicID, document) > 0d;
  }

  @Override
  public String toString() {
    return name;
  }
}
