import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.*;

public class GeneticAlgorithm {
    private static final double A = 10.0;
    private static final double W = 20 * Math.PI;

    private final double[] a;
    private final double[] b;
    private final int[] d;
    private final int N;

    public GeneticAlgorithm(
            double[] a,
            double[] b,
            int[] d,
            int N) {
        this.a = a;
        this.b = b;
        this.d = d;
        this.N = N;
    }

    private double rastrigin(double[] x) {
        double sum = 0.0;
        for (double xi : x) {
            sum += Math.pow(xi, 2) - A * Math.cos(W * xi);
        }
        return A * x.length + sum;
    }

    private int[] computeBinaryLengths() {
        int[] binaryLengths = new int[a.length];

        for (int i = 0; i < a.length; i++) {
            int m = 0;

            while (true) {
                double range = (b[i] - a[i]) * Math.pow(10, d[i]) + 1;

                if (Math.pow(2, m - 1) <= range && range <= Math.pow(2, m)) {
                    break;
                }
                m++;
            }
            binaryLengths[i] = m;
        }

        return binaryLengths;
    }

    private int[] computeNumberOfValues() {
        int[] numberOfValues = new int[a.length];

        for (int i = 0; i < a.length; i++) {
            numberOfValues[i] = (int) ((b[i] - a[i]) * Math.pow(10, d[i]) + 1);
        }

        return numberOfValues;
    }

    private List<Set<String>> createPopulation(int[] binaryLengths, int[] numberOfValues) {
        List<Set<String>> population = new ArrayList<>();
        Random rnd = new Random();

        for (int i = 0; i < binaryLengths.length; i++) {
            int length = binaryLengths[i];
            int maxValue = numberOfValues[i];
            Set<String> individuals = new HashSet<>();

            while (individuals.size() < N) {
                int randomNumber = rnd.nextInt(maxValue);

                String binaryString = String.format("%" + length + "s",
                        Integer.toBinaryString(randomNumber)).replace(' ', '0');

                individuals.add(binaryString);
            }

            population.add(individuals);
        }

        return population;
    }

    private List<String> combineGenesIntoChromosome(List<Set<String>> population) {
        List<String> chromosomes = new ArrayList<>();
        for (int i = 0; i < N; i++) {
            StringBuilder individual = new StringBuilder();

            for (Set<String> gene : population) {
                List<String> list = new ArrayList<>(gene);
                individual.append(list.get(i % list.size()));
            }

            chromosomes.add(individual.toString());
        }
        return chromosomes;
    }

    private List<Double> ratePopulation(List<Set<String>> populations, int[] binaryLengths) {
        List<List<Double>> ratedPopulations = new ArrayList<>();

        for (int index = 0; index < populations.size(); index++) {
            Set<String> population = populations.get(index);
            double a = this.a[index];
            double b = this.b[index];
            int m = binaryLengths[index];

            List<Double> decimalValues = new ArrayList<>();
            for (String individual : population) {
                decimalValues.add((double) Integer.parseInt(individual, 2));
            }

            List<Double> normalizedValues = new ArrayList<>();
            for (Double d : decimalValues) {
                normalizedValues.add(a + ((b - a) * d) / (Math.pow(2, m) - 1));
            }
            ratedPopulations.add(normalizedValues);
        }

        List<List<Double>> transposedPopulations = transpose(ratedPopulations);

        List<Double> evaluationResults = new ArrayList<>();
        transposedPopulations.stream()
                .map(individual -> rastrigin(individual.stream()
                        .mapToDouble(Double::doubleValue)
                        .toArray()))
                .map(result -> new BigDecimal(result)
                        .setScale(2, RoundingMode.HALF_UP)
                        .doubleValue())
                .forEach(evaluationResults::add);


        return evaluationResults;
    }

    public final void optimize() {
        int[] binaryLengths = computeBinaryLengths();
        int[] numberOfValues = computeNumberOfValues();
        List<Set<String>> population = createPopulation(binaryLengths, numberOfValues);
        List<String> chromosomes = combineGenesIntoChromosome(population);
        List<Double> results = ratePopulation(population, binaryLengths);

        for (int i = 0; i < results.size(); i++) {
            System.out.println("Chromosome: " + chromosomes.get(i)
                    + ", Rastrigin value: " + results.get(i));
        }
    }

    private List<List<Double>> transpose(List<List<Double>> matrix) {
        List<List<Double>> transposed = new ArrayList<>();

        int numCols = matrix.get(0).size();

        for (int col = 0; col < numCols; col++) {
            List<Double> newRow = new ArrayList<>();
            for (List<Double> doubles : matrix) {
                newRow.add(doubles.get(col));
            }
            transposed.add(newRow);
        }
        return transposed;
    }

    public static void main(String[] args) {
        double[] a = {-1, -1, -1};
        double[] b = {1, 1, 1};
        int[] d = {1, 1, 1};
        int N = 5;

        GeneticAlgorithm algorithm = new GeneticAlgorithm(a, b, d, N);
        algorithm.optimize();
    }
}