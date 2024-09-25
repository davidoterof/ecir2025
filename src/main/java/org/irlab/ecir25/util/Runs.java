package org.irlab.ecir25.util;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public final class Runs {

  private final List<Map<Integer, List<String>>> runs;

  public Runs(List<Map<Integer, List<String>>> runs) {
    this.runs = runs;
  }

  public static Runs fromListOfPaths(List<String> paths) {
    List<Map<Integer, List<String>>> runs = paths.parallelStream()
                                                 .map(TRECParser::parseRun)
                                                 .sequential()
                                                 .collect(Collectors.toList());
    return new Runs(runs);
  }

  public List<Map<Integer, List<String>>> getRuns() {
    return runs;
  }

  public int size() {
    return runs.size();
  }
}
