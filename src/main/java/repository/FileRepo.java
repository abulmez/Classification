package repository;

import model.Matrix;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
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

    public static Matrix<Double> normalizeData(Matrix<Double> data){
        Matrix<Double> normalizedMatrix = new Matrix<>(data.getRows(),data.getColumns(),0.0);
        ArrayList<Double> columnAverage = new ArrayList<>();
        ArrayList<Double> deviation = new ArrayList<>();
        for(int i=0;i<data.getColumns();i++){
            columnAverage.add(0.0);
            deviation.add(0.0);
        }
        for(int i=0;i<data.getRows();i++)
            for(int j=0;j<data.getColumns();j++)
                columnAverage.set(j,(columnAverage.get(j)+data.get(i,j)));
        for(int i=0;i<data.getColumns();i++){
            columnAverage.set(i,(columnAverage.get(i)/data.getRows()));
        }
        for(int i=0;i<data.getRows();i++)
            for(int j=0;j<data.getColumns();j++)
                deviation.set(j,(deviation.get(j)+Math.pow(data.get(i,j)-columnAverage.get(j),2)));
        for(int i=0;i<data.getColumns();i++){
            deviation.set(i,Math.sqrt(deviation.get(i)/(data.getRows()-1)));
        }
        for(int i=0;i<data.getRows();i++)
            for(int j=0;j<data.getColumns();j++)
                normalizedMatrix.set((data.get(i,j)-columnAverage.get(j))/deviation.get(j),i,j);
        return normalizedMatrix;
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
