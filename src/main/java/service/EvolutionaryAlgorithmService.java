package service;

import DeepCopyUtil.DeepCopy;
import model.Chromosome;
import model.Matrix;
import model.Population;
import repository.FileRepo;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Vector;

public class EvolutionaryAlgorithmService {

    public FileRepo repo;
    private ArrayList<Double> mutationVector;
    public static ArrayList<Double> minHistory;

    public EvolutionaryAlgorithmService(FileRepo repo){
        this.repo = repo;
        minHistory = new ArrayList<>();
    }

    /**
     *  The function initializes the mutation vector with less/higher chance of mutation depending on the mutatianRate parameter
        mutatiaRate - Integer between 0 and 10
     **/
    private void initMutationVector(Integer mutatianRate){
        mutationVector = new ArrayList<>();
        mutationVector.add(0.1*mutatianRate);
        Double remainingMutatianRate = (Double)(1-0.1*mutatianRate);
        mutationVector.add(0.3*remainingMutatianRate);
        mutationVector.add(0.2*remainingMutatianRate);
        mutationVector.add(0.2*remainingMutatianRate);
        mutationVector.add(0.1*remainingMutatianRate);
        mutationVector.add(0.1*remainingMutatianRate);
        mutationVector.add(0.05*remainingMutatianRate);
        mutationVector.add(0.05*remainingMutatianRate);
        cumsum(mutationVector);
    }




    /**
     *  The function determines the value of a given probability using a provided simple random variable
     * @param list  The distribution of the random variable
     * @param random Double (probability between 0 and 1)
     * @return  Returns the corresponding value for the given probability
     */
    private Integer vasInv(List<Double> list, Double random){
        for(int i=0;i<list.size();i++){
            if(random<list.get(i)){
                return i;
            }
        }
        return -1;
    }

    /**
     * The function determines the value associated with a chromosome
     * @param c Chromosome
     * @return  Returns the values of the provided Chromosome
     */
    private double fitness(Matrix<Double> data,Matrix<Double> results,Chromosome c, Integer resultColumn){
        Double totalError = 0.0;
        for(int i=0;i<data.getRows();i++) {
            Double determinedValue = c.getValue(0);
            for (int j = 0; j < data.getColumns(); j++) {
                determinedValue += data.get(i,j)*c.getValue(j+1);
            }
            totalError += Math.abs(results.get(i,resultColumn)-sigmoidFunction(determinedValue));
        }
        return totalError;
    }

    /**
     * The function determines the cumulative sum of a list of values
     * @param list List of Double values
     */
    private void cumsum(List<Double> list){
        for(int j=1;j<list.size();j++){
            list.set(j,list.get(j)+list.get(j-1));
        }
    }

    private Double sigmoidFunction(Double x){
        return 1 / (1 + Math.exp(-x));
    }

    /**
     * The function determines the minimum value in a population
     * @param p Population
     * @return Returns the minimum found value
     */
    private ArrayList<Double> evaluatePopulation(Matrix<Double> data,Matrix<Double> results,Population p, Integer resultColumn){
        ArrayList<Double> errorList = new ArrayList<>();
        for(Chromosome c:p.getChromosomesList()){
            errorList.add(fitness(data,results,c,resultColumn));
        }
        return errorList;
    }




    /**
     * The function mutates a Chromosome
     * @param offspring The Chromosome to be mutated
     */
    private  void mutate(Chromosome offspring){
        for(int i=0;i<offspring.getValues().size();i++){
            Double mutationOffset =  (double) (vasInv(mutationVector, Math.random()));
            mutationOffset/=1000;
            Boolean mutationDirection = Math.random()<0.5;
            Double mutatedValue;
            if(mutationDirection)
                mutatedValue =  offspring.getValue(i)+mutationOffset;
            else
                mutatedValue = offspring.getValue(i)-mutationOffset;
            offspring.setValue(i,mutatedValue);
        }
    }

    /**
     * The function creates a new Chromosome according to two other Chromosomes
     * @param dad  One of the parent Chromosome
     * @param mom  The other parent Chromosome
     * @return  Return the offspring of the 2 provided chromosomes
     */
    private Chromosome offspring(Chromosome dad, Chromosome mom){
        Vector<Double> newValues = new Vector<>();
        for(int i=0;i<dad.getValues().size();i++){
            newValues.add((dad.getValue(i)+mom.getValue(i))/2);
        }
        return new Chromosome(newValues.toArray(new Double[0]));
    }

    private Integer getBestChromosomeInPopulation(Population p, ArrayList<Double> errorList){
        Double min= Double.MAX_VALUE;
        Integer minPosition = 0;
        for(int i=0;i<errorList.size();i++){
            if(errorList.get(i)<min){
                minPosition = i;
                min = errorList.get(i);
            }
        }
        return minPosition;
    }

    private ArrayList<Double> sortChances(ArrayList<Double> errorList,ArrayList<Double> chance){
        for(int i=0;i<chance.size();i++){
            for(int j=0;j<chance.size();j++){
                if(i!=j && errorList.get(i)<errorList.get(j) && chance.get(i)<chance.get(j)){
                    Double aux = chance.get(i);
                    chance.set(i,chance.get(j));
                    chance.set(j,aux);
                }
            }
        }
        return chance;
    }

    /**
     * The function determines the minimum value in an 1 or 2 dimensional array using an evolutionary algorithm
     * @param numberOfGenerations
     * @param populationSize    The size of every population in each generation
     * @param mutatianRate     The mutation rate of the population
     * @return  Return the found minimum
     */
    public ArrayList<Double> solve(Integer numberOfGenerations, Integer populationSize, Integer mutatianRate, Integer resultColumn){
        initMutationVector(10-mutatianRate);
        Matrix<Double> normalizedData = FileRepo.normalizeData(repo.getDataMatrix());
        Matrix<Double> normalizedResults = repo.getResultsMatrix();
        //Matrix<Double> normalizedResults = normalizeData(repo.getResultsMatrix());
        Population p = new Population();
        for(int i=0;i<populationSize;i++) {
            ArrayList<Double> chromosomeValues = new ArrayList<>();
            for(int j=0;j<normalizedData.getColumns()+1;j++){
                Double d = Math.random();
                chromosomeValues.add(d);
            }
            p.addChromosome(new Chromosome(chromosomeValues.toArray(new Double[0])));
        }
        for(int i=0;i<numberOfGenerations;i++){
            Population offsprings = new Population();
            ArrayList<Double> chance = new ArrayList<>();
            ArrayList<Double> errorList = evaluatePopulation(normalizedData,normalizedResults,p,resultColumn);
            Double errorSum = errorList.stream().mapToDouble(Double::doubleValue).sum();
            for(Double elem:errorList){
                chance.add(elem/errorSum);
            }
            sortChances(errorList,chance);
            cumsum(chance);
            for(int j=0;j<populationSize-1;j++) {
                Chromosome dad = p.getChromosome(vasInv(chance, Math.random()));
                Chromosome mom = p.getChromosome(vasInv(chance, Math.random()));
                Chromosome offspring = offspring(dad,mom);
                mutate(offspring);
                offsprings.addChromosome(offspring);
            }
            Integer bestChromosomePosition = getBestChromosomeInPopulation(p,errorList);
            Chromosome bestFromLastPopulation = (Chromosome) DeepCopy.copy(p.getChromosome(bestChromosomePosition));
            offsprings.addChromosome(bestFromLastPopulation);
            p = offsprings;
        }
        ArrayList<Double> errorList = evaluatePopulation(normalizedData,normalizedResults,p,resultColumn);
        Integer bestChromosomePosition = getBestChromosomeInPopulation(p,errorList);
        //System.out.println("Eroarea totala pentru coeficientii determinati: "+errorList.get(bestChromosomePosition));
        return p.getChromosome(bestChromosomePosition).getValues();
    }

    public Double testAccuracy(ArrayList<Double> factor, Integer resultColumn){
        Integer ok = 0;
        Matrix<Double> normalizedTestData = FileRepo.normalizeData(repo.getTestDataMatrix());
        Matrix<Double> normalizedTestDataResults = repo.getTestResultsMatrix();
        for(int i=0;i<normalizedTestData.getRows();i++){
            Double resultValue = factor.get(0);
            for(int j=1;j<factor.size();j++){
                resultValue += factor.get(j)*normalizedTestData.get(i,j-1);
            }

            resultValue = sigmoidFunction(resultValue);
            if(resultValue<=(1.0/3))
                resultValue = 0.0;
            else if(resultValue>(1.0/3) && resultValue<=(2.0/3) )
                resultValue = 0.5;
            else
                resultValue = 1.0;
            if(resultValue.equals(normalizedTestDataResults.get(i,resultColumn)))
                ok += 1;
        }
        return ok/(double)normalizedTestData.getRows();
    }

    /**
     * The function determines the minimum value in an 1 or 2 dimensional array by browsing the array
     * @return Return the minimum
     */


}
