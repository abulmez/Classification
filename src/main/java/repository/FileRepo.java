package repository;

import model.Matrix;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

public class FileRepo {

    private Matrix<Double> dataMatrix = new Matrix<>(1000,6,0.0);
    private Matrix<Double> resultsMatrix = new Matrix<>(1000, 1,0.0);
    private Matrix<Double> testDataMatrix = new Matrix<>(1000,6,0.0);
    private Matrix<Double> testResultsMatrix = new Matrix<>(1000, 1,0.0);

    public FileRepo(String trainingDataFileName, String testingDataFileName){
        readDataFromFile(dataMatrix,resultsMatrix,trainingDataFileName);
        readDataFromFile(testDataMatrix,testResultsMatrix,testingDataFileName);
    }

    private void readDataFromFile(Matrix<Double> data, Matrix<Double> result, String fileName){
        try{
            Integer counter = 0;
            ClassLoader classLoader = getClass().getClassLoader();
            File file =
                    new File(classLoader.getResource(fileName).getFile());
            Scanner sc = new Scanner(file);
            String line;
            try {
                while ((line = sc.nextLine()) != null) {

                    String[] values = line.split(",");
                    for (int j = 0; j < 6; j++) {
                        data.set(Double.parseDouble(values[j]), counter, j);

                    }
                    result.set(Double.parseDouble(values[6]), counter, 0);
                    counter += 1;
                }
            }
            catch (Exception e){
                data.setRows(counter);
                result.setRows(counter);
            }
        }
        catch (FileNotFoundException | NullPointerException e) {
            e.printStackTrace();
        }
    }

    public Matrix<Double> getDataMatrix() {
        return dataMatrix;
    }

    public Matrix<Double> getResultsMatrix() {
        return resultsMatrix;
    }

    public Matrix<Double> getTestDataMatrix() {
        return testDataMatrix;
    }

    public Matrix<Double> getTestResultsMatrix() {
        return testResultsMatrix;
    }
}
