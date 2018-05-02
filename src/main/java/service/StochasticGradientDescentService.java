package service;

import model.Matrix;
import repository.FileRepo;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

public class StochasticGradientDescentService {

    FileRepo repo;

    public StochasticGradientDescentService(FileRepo repo) {
        this.repo = repo;
    }

    public ArrayList<Double> getError(Matrix<Double> data,Matrix<Double> results,ArrayList<Double> factor, Integer resultColumn){

        ArrayList<Double> errorArray = new ArrayList<>();
        for(int i=0;i<data.getRows();i++) {
            Double currentValue = factor.get(0);
            for (int j = 1; j < factor.size(); j++) {
                currentValue += data.get(i, j - 1) * factor.get(j);
            }
            errorArray.add(results.get(i,resultColumn)-sigmoidFunction(currentValue));
        }
        return errorArray;
    }

    private Double sigmoidFunction(Double x){
        return 1 / (1 + Math.exp(-x));
    }

    public ArrayList<Double> solve(Integer numberOfIterations, Double learningRate, Integer resultColumn){

        Matrix<Double> normalizedData = FileRepo.normalizeData(repo.getDataMatrix());
        Matrix<Double> normalizedResults = repo.getResultsMatrix();
        learningRate/=normalizedData.getRows();
        ArrayList<Integer> indexes = new ArrayList<>();
        ArrayList<Double> factor = new ArrayList<>();
        for(int i=0;i<normalizedData.getRows();i++){
            indexes.add(i);
        }
        Collections.shuffle(indexes);
        for(int i=0;i<normalizedData.getColumns()+1;i++){
            factor.add(Math.random());
        }
        ArrayList<Double> errors;
        for(int i=0;i<numberOfIterations;i++) {
            errors = getError(normalizedData,normalizedResults,factor, resultColumn);
            for (int j = 1; j < factor.size(); j++) {
                Double gradient = 0.0;
                for(int k=0;k<errors.size();k++){
                    gradient+=normalizedData.get(k,j-1)*errors.get(k);
                }
                factor.set(j, (factor.get(j) - learningRate*gradient));
            }
            factor.set(0, (factor.get(0) - learningRate*(errors.stream().mapToDouble(Double::doubleValue).sum())));
        }
        return factor;
    }
}
