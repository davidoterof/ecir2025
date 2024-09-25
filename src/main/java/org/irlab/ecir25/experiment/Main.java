package org.irlab.ecir25.experiment;

import org.apache.commons.cli.*;

public class Main {

  public static void main(String[] args) {
    CommandLine parsedArgs = parse(args);

    if (parsedArgs == null) {
      return;
    }

    if (parsedArgs.getOptionValue("e").equals("simulation")) {
      ExperimentCompareTests.run(parsedArgs);
    } else if (parsedArgs.getOptionValue("e").equals("mq")) {
      ExperimentMillionQuery.run(parsedArgs);
    } else {
      System.out.println("Wrong experiment type!");
    }
  }

  private static Options getOptions() {
    Options options = new Options();

    Option experiment = Option.builder("e")
                              .longOpt("experiment")
                              .desc("experiment to perform: simulation or mq")
                              .required()
                              .hasArgs()
                              .argName("type")
                              .build();

    Option collection = Option.builder("c")
                              .longOpt("collection")
                              .desc("name of the collection")
                              .required(false)
                              .hasArgs()
                              .argName("name")
                              .build();

    Option output = Option.builder("o")
                          .longOpt("output")
                          .desc("path to write experimental results")
                          .required(true)
                          .hasArgs()
                          .argName("path")
                          .build();

    Option runs = Option.builder("r")
                        .longOpt("runs")
                        .desc("path to runs location")
                        .required(true)
                        .hasArgs()
                        .argName("path")
                        .build();

    Option qrels = Option.builder("q")
                         .longOpt("qrels")
                         .desc("path to qrels file")
                         .required(true)
                         .hasArgs()
                         .argName("path")
                         .build();

    options.addOption(experiment);
    options.addOption(collection);
    options.addOption(output);
    options.addOption(qrels);
    options.addOption(runs);
    options.addOption(runs);

    return options;
  }

  private static CommandLine parse(String[] args) {
    Options options = getOptions();

    CommandLineParser parser = new DefaultParser();
    try {
      return parser.parse(getOptions(), args);
    } catch (ParseException e) {
      HelpFormatter formatter = new HelpFormatter();
      formatter.printHelp("ecir25", options, true);
      return null;
    }
  }
}
