This repo holds the code to reproduce the experiments of the following paper:

#### Towards Reliable Testing for Multiple Information Retrieval System Comparisons @ ECIR '25 ([Link](https://arxiv.org/abs/2501.03930))

To run these experiments, you will need maven and the JDK 11 installed.
Then, you should compile the project with

```
mvn clean package
```

This will create a fatjar under the ```target/``` folder. Execute this jar to launch the experiments.

## Simulation experiments

To run the simulation experiments, use the option ```-e simulation```. Use
```java -jar target/ecir25-1.0-jar-with-dependencies.jar -h```
to see the rest of the options.

It will create several directories under the output path, each directory
for a different combination of the parameters of the experiment (number of systems,
sample size, metric, etc.). In each directory, there will be a file for each
run in the collection. These files then will be use to plot the FWER and the
power of each test (more on how to plot the results in the last section of this
readme)

BE CAREFUL when running these experiments, as it will take every core available on the machine.

You can tune the number of simulations performed so that the experiments do not take so long.
This parameter is in the file ```ExperimentCompareTests.java```.

## Million Query experiment

To run the Million Query experiment, use the option ```-e mq```. You will need
to obtain the runs and the qrels of this dataset.

This experiment will create a lot of different files in the ouput folder, so use
a different one than in the previous experiment.

BE CAREFUL when running these experiments, as it will take every core available on the machine.

You can tune the static variables in ```ExperimentMillionquery``` so that the experiments
do not take so long.

## Plotting

We also release scripts to plot same figures as the paper. To run them, the following dependencies are needed:

```
pandas
numpy
matplotlib
```

All scripts are under ```src/main/py```.
