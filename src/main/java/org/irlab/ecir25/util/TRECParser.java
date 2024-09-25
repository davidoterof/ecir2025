package org.irlab.ecir25.util;

import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.objects.Object2DoubleLinkedOpenHashMap;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public final class TRECParser {

  private static final Logger LOG = LogManager.getLogger();

  private TRECParser() {
    throw new AssertionError("Must not instantiate an element of this class");
  }

  public static Map<Integer, Map<String, Double>> parseQrels(final String file) {
    final Path qrelsPath = Paths.get(file);
    final Map<Integer, Map<String, Double>> qrels = new Int2ObjectOpenHashMap<>();

    if (Files.notExists(qrelsPath)) {
      LOG.error("qrels file does not exist: {}", file);
      return qrels;
    }

    try (BufferedReader reader = Files.newBufferedReader(qrelsPath)) {

      String line;
      while ((line = reader.readLine()) != null) {
        final String[] splits = line.split("\\s+");

        final int queryID = Integer.parseInt(splits[0]);
        final String docID = splits[2];
        final double relevance = Double.parseDouble(splits[3]);

        if (!qrels.containsKey(queryID)) {
          // Use a linked map to maintain the order of the qrels as are read from disk.
          qrels.put(queryID, new Object2DoubleLinkedOpenHashMap<>());
        }

        qrels.get(queryID).put(docID, relevance);
      }
    } catch (final IOException e) {
      LOG.error(e.getMessage());
    }

    return qrels;
  }

  public static Map<Integer, List<String>> parseRun(String file) {
    final Map<Integer, List<String>> finalRun = new Int2ObjectOpenHashMap<>();

    try (BufferedReader reader = Files.newBufferedReader(Paths.get(file))) {
      String line;
      while ((line = reader.readLine()) != null) {
        final String[] parts = line.trim().split("\\s+");
        final int queryID = Integer.parseInt(parts[0]);

        if (!finalRun.containsKey(queryID)) {
          finalRun.put(queryID, new ArrayList<>());
        }

        final List<String> queryRankingDocs = finalRun.get(queryID);

        final String doc = parts[2];
        queryRankingDocs.add(doc);
      }

    } catch (final IOException e1) {
      LOG.error(e1.getMessage());
    }

    return finalRun;
  }
}
