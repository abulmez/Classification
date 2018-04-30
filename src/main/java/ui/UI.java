package ui;

import model.Matrix;
import repository.FileRepo;
import service.EvolutionaryAlgorithmService;
import java.util.ArrayList;
import java.util.Scanner;

public class UI {

    private EvolutionaryAlgorithmService evolutionaryAlgorithmService;

    private Scanner sc;
    public UI(){
        sc = new Scanner(System.in);
    }

    /**
     * Displays the main menu of the application
     */
    public void displayMainMenu() {

        int x=2;
        FileRepo repo = new FileRepo("column_3C_weka_data.arff", "column_3C_weka_test.arff");
        evolutionaryAlgorithmService = new EvolutionaryAlgorithmService(repo);



        Boolean merge = true;
        while(merge) {

            System.out.println("Alegeti algoritmul de antrenare:");
            System.out.println("1.Gradient descendent");
            System.out.println("2.Algoritm evolutiv");
            System.out.println("3.Iesire");
            x = sc.nextInt();


            switch (x) {


                case 2: {
                    // System.out.print("done");
                    System.out.println("Dati numarul de generatii:");
                    int generationNumber = sc.nextInt();
                    System.out.println("Dati marimea populatiei:");
                    int populationSize = sc.nextInt();
                    System.out.println("Dati factorul de mutatie: (valoare intre 0 si 10)");
                    int mutatianRate = sc.nextInt();
                    int resultColumn = 0;
                    ArrayList<Double> result = evolutionaryAlgorithmService.solve(generationNumber, populationSize, mutatianRate, resultColumn );
                    System.out.println("Coeficientii determinati:");
                    for (Double rez : result) {
                        System.out.print(rez + " ");
                    }
                    System.out.println();
                    System.out.println("Performanta algoritmului: " + (int)(evolutionaryAlgorithmService.testAccuracy(result,resultColumn)*100)+"%");
                    break;
                }


                case 1: {
                    /*
                    System.out.println("Dati rezultatul pentru care se va efectua regresia:");
                    System.out.println("1.motor_UPDRS");
                    System.out.println("2.total_UPDRS");
                    int resultColumn = sc.nextInt();
                    System.out.println("Dati gradul de invatare:");
                    Double learningRate = sc.nextDouble();
                    System.out.println("Dati numarul de iteratii:");
                    int numberOfIterations = sc.nextInt();
                    ArrayList<Double> result = stochasticGradientDescentService.solve(numberOfIterations, learningRate, resultColumn-1);
                    System.out.println("Coeficientii determinati:");
                    for (Double rez : result) {
                        System.out.print(rez + " ");
                    }
                    System.out.println();
                    System.out.println("Valoarea reala | Valoarea Determinata | Eroare");
                    System.out.println();
                    Matrix<Double> evaluationResult = evolutionaryAlgorithmService.test(result, resultColumn - 1);
                    Double error = 0.0;
                    for (int i = 0; i < evolutionaryAlgorithmService.repo.getTestResultMatrixNumberOfRows(); i++) {
                        System.out.println(evaluationResult.get(i, 0) + " " + evaluationResult.get(i, 1) + " " + Math.abs(evaluationResult.get(i, 0) - evaluationResult.get(i, 1)));
                        error += Math.abs(evaluationResult.get(i, 0) - evaluationResult.get(i, 1));
                    }
                    System.out.println();
                    System.out.println("Eroarea medie: "+error/evolutionaryAlgorithmService.repo.getTestResultMatrixNumberOfRows());
                    */
                    break;
                }
                case 3: {
                    merge = false;
                    break;
                }

                default:
                    System.out.println("Optiune invalida! Mai incercati!");
            }
        }


    }
}
